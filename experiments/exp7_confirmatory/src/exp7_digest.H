// Experiment 7 - streaming digests over returned floating-point output.
//
// Strict C++11, no dependency outside the standard library.
//
// Why this file exists.  The candidate's random stream is a function of the engine alone,
// so two builds of this probe against different standard libraries must produce *identical
// bytes* for every CANDIDATE configuration.  Testing that with summary statistics and
// interval overlap would be a weaker claim than the one being made.  Each configuration
// therefore carries a digest over the IEEE-754 bit patterns of every returned component,
// in order, followed by a canonical serialization of every counter; the analysis asserts
// exact equality of those digests across standard libraries rather than agreement of rates.
//
// Two digests are emitted for each stream.  SHA-256 is the one the claim rests on.  FNV-1a
// is carried alongside because it is cheap, is computed over the same byte sequence, and
// gives an independent signal if a build's SHA-256 were itself miscompiled -- the two would
// then disagree about *which* streams match.
//
// The bit patterns, not the values, are digested: a NaN payload, a signed zero and an
// infinity are all distinguishable this way, and all three occur in this experiment.

#ifndef EXP7_DIGEST_H
#define EXP7_DIGEST_H

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace exp7
{

// ---------------------------------------------------------------------------
// FNV-1a, 64 bit (Fowler-Noll-Vo).  Offset basis and prime as published.
// ---------------------------------------------------------------------------
class Fnv1a64
{
  public:
    Fnv1a64() : m_h(14695981039346656037ULL) {}

    void update(const void *data, size_t n)
    {
        const unsigned char *p= static_cast<const unsigned char *>(data);
        for (size_t i= 0; i < n; ++i)
        {
            m_h^= static_cast<uint64_t>(p[i]);
            m_h*= 1099511628211ULL;
        }
    }

    uint64_t value() const { return m_h; }

    std::string hex() const
    {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(m_h));
        return std::string(buf);
    }

  private:
    uint64_t m_h;
};

// ---------------------------------------------------------------------------
// SHA-256 (FIPS 180-4).  Written out here rather than taken from a library so that the
// probe keeps its "no dependency outside the standard library" property; `make selftest`
// runs the published known-answer vectors against it before any production phase.
// ---------------------------------------------------------------------------
class Sha256
{
  public:
    Sha256() { reset(); }

    void reset()
    {
        m_len= 0;
        m_fill= 0;
        m_h[0]= 0x6a09e667u; m_h[1]= 0xbb67ae85u; m_h[2]= 0x3c6ef372u; m_h[3]= 0xa54ff53au;
        m_h[4]= 0x510e527fu; m_h[5]= 0x9b05688cu; m_h[6]= 0x1f83d9abu; m_h[7]= 0x5be0cd19u;
    }

    void update(const void *data, size_t n)
    {
        const unsigned char *p= static_cast<const unsigned char *>(data);
        m_len+= static_cast<uint64_t>(n);
        while (n > 0)
        {
            const size_t take= (64 - m_fill < n) ? (64 - m_fill) : n;
            std::memcpy(m_block + m_fill, p, take);
            m_fill+= take;
            p+= take;
            n-= take;
            if (m_fill == 64)
            {
                transform(m_block);
                m_fill= 0;
            }
        }
    }

    /// Finalizes a copy, so the object may keep being fed afterwards.
    std::string hex() const
    {
        Sha256 c(*this);
        unsigned char out[32];
        c.finish(out);
        static const char *kHex= "0123456789abcdef";
        std::string s;
        s.reserve(64);
        for (int i= 0; i < 32; ++i)
        {
            s.push_back(kHex[(out[i] >> 4) & 0xF]);
            s.push_back(kHex[out[i] & 0xF]);
        }
        return s;
    }

    /// SHA-256 of a whole file, or an empty string when it cannot be read.  Used for the
    /// run-time manifest, which records the hash of each output file as written rather
    /// than the hash of whatever happens to be on disk when the analysis runs.
    static std::string file(const std::string &path)
    {
        FILE *f= std::fopen(path.c_str(), "rb");
        if (!f)
            return std::string();
        Sha256 h;
        std::vector<unsigned char> buf(1 << 16);
        for (;;)
        {
            const size_t got= std::fread(&buf[0], 1, buf.size(), f);
            if (got > 0)
                h.update(&buf[0], got);
            if (got < buf.size())
                break;
        }
        std::fclose(f);
        return h.hex();
    }

  private:
    static uint32_t ror(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

    void finish(unsigned char out[32])
    {
        const uint64_t bits= m_len * 8ULL;
        unsigned char pad= 0x80;
        update(&pad, 1);
        pad= 0x00;
        while (m_fill != 56)
            update(&pad, 1);
        unsigned char lenbe[8];
        for (int i= 0; i < 8; ++i)
            lenbe[i]= static_cast<unsigned char>((bits >> (56 - 8 * i)) & 0xFF);
        update(lenbe, 8);
        for (int i= 0; i < 8; ++i)
        {
            out[4 * i + 0]= static_cast<unsigned char>((m_h[i] >> 24) & 0xFF);
            out[4 * i + 1]= static_cast<unsigned char>((m_h[i] >> 16) & 0xFF);
            out[4 * i + 2]= static_cast<unsigned char>((m_h[i] >> 8) & 0xFF);
            out[4 * i + 3]= static_cast<unsigned char>(m_h[i] & 0xFF);
        }
    }

    void transform(const unsigned char block[64])
    {
        static const uint32_t K[64]= {
            0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
            0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
            0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
            0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
            0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
            0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
            0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
            0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
            0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
            0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
            0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};

        uint32_t w[64];
        for (int i= 0; i < 16; ++i)
            w[i]= (static_cast<uint32_t>(block[4 * i]) << 24) |
                  (static_cast<uint32_t>(block[4 * i + 1]) << 16) |
                  (static_cast<uint32_t>(block[4 * i + 2]) << 8) |
                  static_cast<uint32_t>(block[4 * i + 3]);
        for (int i= 16; i < 64; ++i)
        {
            const uint32_t s0= ror(w[i - 15], 7) ^ ror(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const uint32_t s1= ror(w[i - 2], 17) ^ ror(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i]= w[i - 16] + s0 + w[i - 7] + s1;
        }

        uint32_t a= m_h[0], b= m_h[1], c= m_h[2], d= m_h[3];
        uint32_t e= m_h[4], f= m_h[5], g= m_h[6], h= m_h[7];
        for (int i= 0; i < 64; ++i)
        {
            const uint32_t S1= ror(e, 6) ^ ror(e, 11) ^ ror(e, 25);
            const uint32_t ch= (e & f) ^ ((~e) & g);
            const uint32_t t1= h + S1 + ch + K[i] + w[i];
            const uint32_t S0= ror(a, 2) ^ ror(a, 13) ^ ror(a, 22);
            const uint32_t mj= (a & b) ^ (a & c) ^ (b & c);
            const uint32_t t2= S0 + mj;
            h= g; g= f; f= e; e= d + t1; d= c; c= b; b= a; a= t1 + t2;
        }
        m_h[0]+= a; m_h[1]+= b; m_h[2]+= c; m_h[3]+= d;
        m_h[4]+= e; m_h[5]+= f; m_h[6]+= g; m_h[7]+= h;
    }

    uint32_t m_h[8];
    unsigned char m_block[64];
    size_t m_fill;
    uint64_t m_len;
};

// ---------------------------------------------------------------------------
// A stream digest that feeds both hashes the same byte sequence.
//
// Floating-point values enter as their raw object representation, so a signed zero, an
// infinity and a NaN payload are all distinguishable.  Integers enter as fixed-width
// little-endian, so the digest does not depend on the host's integer width choices.
// ---------------------------------------------------------------------------
class StreamDigest
{
  public:
    StreamDigest() : m_bytes(0) {}

    void bytes(const void *p, size_t n)
    {
        m_sha.update(p, n);
        m_fnv.update(p, n);
        m_bytes+= static_cast<uint64_t>(n);
    }

    /// Raw object representation of one float or double.
    template <typename T> void real(T v)
    {
        unsigned char b[sizeof(T)];
        std::memcpy(b, &v, sizeof(T));
        bytes(b, sizeof(T));
    }

    /// Little-endian 64-bit integer, whatever the host's byte order.
    void u64(uint64_t v)
    {
        unsigned char b[8];
        for (int i= 0; i < 8; ++i)
            b[i]= static_cast<unsigned char>((v >> (8 * i)) & 0xFF);
        bytes(b, 8);
    }

    void i64(long long v) { u64(static_cast<uint64_t>(v)); }

    /// A name, so that the same integer under two different counters cannot collide.
    void label(const char *s) { bytes(s, std::strlen(s) + 1); }

    void counter(const char *name, long long v)
    {
        label(name);
        i64(v);
    }

    uint64_t byteCount() const { return m_bytes; }
    std::string sha256() const { return m_sha.hex(); }
    std::string fnv1a() const { return m_fnv.hex(); }

  private:
    Sha256 m_sha;
    Fnv1a64 m_fnv;
    uint64_t m_bytes;
};

} // namespace exp7

#endif // EXP7_DIGEST_H
