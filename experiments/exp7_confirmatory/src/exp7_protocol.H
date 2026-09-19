// GENERATED FILE -- DO NOT EDIT.
//
// Written by src/gen_protocol_header.py from config/protocol.json.  It carries the
// parts of the frozen protocol that the probe itself has to act on: the quantile
// levels and their order-statistic index brackets, the upper-tail thresholds, the
// incomplete-beta switch point, the frozen sample sizes and the matrix.  Everything
// else the protocol fixes is a decision rule, and decision rules belong to the
// analysis, which reads the JSON directly.
//
// config/protocol.json SHA-256: e336f6a392bdf2f370a4f368f5816a1847b0b6d41150c9c3b4ff5a11cf608971

#ifndef EXP7_PROTOCOL_H
#define EXP7_PROTOCOL_H

namespace exp7
{
namespace protocol
{

static const char *const kProtocolVersion= "1.2.0";
static const char *const kProtocolSha256= "e336f6a392bdf2f370a4f368f5816a1847b0b6d41150c9c3b4ff5a11cf608971";

/// Intended attempts per seed, by phase, exactly as PROTOCOL.md Sec. 4 freezes them.
static const long long kNScalar= 1000000;
static const long long kNMechanism= 1000000;
static const long long kNConditioning= 1000000;
static const long long kNLoader= 100000;

/// Quantile levels of F2/F3, and the exact 0-based order-statistic bracket the
/// protocol froze for each of them at n = kNScalar.  The probe emits the sample
/// values at these indices; the analysis owns the coverage rule.
static const int kNumQuantileLevels= 5;
static const double kQuantileLevels[kNumQuantileLevels]= {
    0.5, 0.9, 0.99, 0.999, 0.9999};
static const long long kQuantileLoIndex[kNumQuantileLevels]= {
    498711, 899226, 989742, 998917, 999872};
static const long long kQuantileHiIndex[kNumQuantileLevels]= {
    501288, 900772, 990255, 999080, 999925};
static const double kQuantileAchievedCoverage[kNumQuantileLevels]= {
    0.9900338398639903, 0.9900243176826826, 0.9900610423083096, 0.9900855825109875, 0.9920353224105682};

/// F4 upper-tail thresholds: exceedance counts of Z above z0 = -log(q0).
static const int kNumTailQ0= 3;
static const double kTailQ0[kNumTailQ0]= {
    0.01, 0.001, 0.0001};

/// Switch point of the two-piece log I_W(a, 3/2) evaluator, in log W.
static const double kLogQSwitch= -30.0;

/// Amendment 1.1.0: a returned draw whose radius differs from the exact value by
/// more than this fraction is FINITE_BUT_WRONG -- a loss, not a success.  The
/// threshold is 2^-(digits/2), i.e. more than half the significand of its type gone.
static const double kMaxRelErrorDouble= 1.0536712127723509e-08;
static const double kMaxRelErrorFloat= 0.000244140625;

/// Amendment 1.2.0: the audit stratum sampling rates and the decision-relevance
/// margin are protocol data, not implementation details.  A rate of 1 means the
/// stratum is adjudicated in full.
static const double kAuditMarginLogUnits= 20.0;
static const int kNumAuditRates= 7;
static const char *const kAuditStratumNames[kNumAuditRates]= {
    "within_margin_of_a_type_limit",
    "method_disagreement",
    "finite_but_wrong_candidate",
    "avoidable_loss_candidate",
    "subnormal_or_zero_denominator_representable_target",
    "unambiguous_failure_beyond_margin",
    "uniform_sample"};
static const double kAuditStratumRates[kNumAuditRates]= {
    1.0, 1.0, 1.0, 1.0, 1.0, 0.001, 0.0001};

/// The matrix, restated here so that the hand-written declarations in
/// exp7_common.H can be checked against the frozen JSON at start-up.
static const int kNumKappa= 13;
static const double kKappaLadder[kNumKappa]= {
    0.5001, 0.501, 0.505, 0.51, 0.55, 0.6, 0.75, 1.0, 1.25, 1.49, 1.5, 2.0, 5.0};
static const int kNumProductionSeeds= 5;
static const unsigned kProductionSeeds[kNumProductionSeeds]= {
    7001u, 7002u, 7003u, 7004u, 7005u};
static const int kNumPerformanceSeeds= 5;
static const unsigned kPerformanceSeeds[kNumPerformanceSeeds]= {
    7006u, 7007u, 7008u, 7009u, 7010u};

} // namespace protocol
} // namespace exp7

#endif // EXP7_PROTOCOL_H
