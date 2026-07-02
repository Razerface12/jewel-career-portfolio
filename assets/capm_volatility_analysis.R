# =============================================================================
#  CAPM BETA & VOLATILITY ANALYSIS  -  A.G. Barr (BAG.L) vs FTSE 250
#  Author: Jewel Thomas  |  Portfolio project (equity risk / quantitative)
# -----------------------------------------------------------------------------
#  WHAT THIS DOES
#   1. Loads ~3 years of weekly prices for a stock and a market index.
#   2. Converts prices to log returns.
#   3. Measures risk: annualised volatility + a 26-week rolling volatility.
#   4. Estimates the CAPM: regresses the stock's EXCESS return on the market's
#      EXCESS return to obtain beta (systematic risk), alpha and R-squared.
#   5. Saves two charts and prints a clean results summary.
#
#  HOW TO RUN
#   * LIVE DATA (default): needs internet + the 'quantmod' package. Pulls prices
#     straight from Yahoo Finance, so results reflect the latest market data.
#   * OFFLINE / REPRODUCIBLE: set USE_LIVE_DATA <- FALSE to read the bundled
#     'sample_prices.csv' instead (base R only, runs anywhere).
#
#  NOTE ON THE SAMPLE FILE: sample_prices.csv is an ILLUSTRATIVE dataset so the
#  script is reproducible without a data feed. Run with USE_LIVE_DATA <- TRUE to
#  reproduce the analysis on real BAG.L / FTSE 250 prices.
# =============================================================================

USE_LIVE_DATA <- TRUE          # <- set FALSE to use the bundled sample CSV
RISK_FREE_ANNUAL <- 0.04       # annual risk-free rate (UK ~ base rate proxy)
ROLL_WINDOW <- 26              # rolling volatility window, in weeks
PERIODS_PER_YEAR <- 52         # weekly data -> 52 periods per year

# -----------------------------------------------------------------------------
# 1. LOAD DATA  ->  data.frame with columns: Date, stock, market
# -----------------------------------------------------------------------------
load_prices <- function() {
  if (USE_LIVE_DATA) {
    # requires: install.packages("quantmod")
    if (!requireNamespace("quantmod", quietly = TRUE))
      stop("Package 'quantmod' not installed. install.packages('quantmod') or set USE_LIVE_DATA <- FALSE")
    library(quantmod)
    # BAG.L = A.G. Barr; ^FTMC = FTSE 250 index (change tickers to re-use)
    getSymbols(c("BAG.L", "^FTMC"), src = "yahoo", from = Sys.Date() - 365*3, periodicity = "weekly")
    stock  <- Cl(get("BAG.L"))          # weekly closing price
    market <- Cl(get("FTMC"))
    merged <- na.omit(merge(stock, market))
    data.frame(Date = index(merged),
               stock = as.numeric(merged[, 1]),
               market = as.numeric(merged[, 2]))
  } else {
    df <- read.csv("sample_prices.csv", stringsAsFactors = FALSE)
    data.frame(Date = as.Date(df$Date), stock = df[, 2], market = df[, 3])
  }
}

px <- load_prices()

# -----------------------------------------------------------------------------
# 2. LOG RETURNS  (log returns are additive over time -> cleaner for stats)
# -----------------------------------------------------------------------------
log_ret <- function(p) diff(log(p))
rs <- log_ret(px$stock)      # stock weekly log returns
rm <- log_ret(px$market)     # market weekly log returns
rf <- RISK_FREE_ANNUAL / PERIODS_PER_YEAR   # per-period risk-free rate

# -----------------------------------------------------------------------------
# 3. VOLATILITY  (annualised = weekly sd * sqrt(52))
# -----------------------------------------------------------------------------
ann_vol <- function(r) sd(r) * sqrt(PERIODS_PER_YEAR)
vol_stock  <- ann_vol(rs)
vol_market <- ann_vol(rm)

# rolling annualised volatility (risk through time, not just an average)
rolling_vol <- function(r, w) {
  out <- rep(NA_real_, length(r))
  for (i in w:length(r)) out[i] <- sd(r[(i - w + 1):i]) * sqrt(PERIODS_PER_YEAR)
  out
}
rv_stock  <- rolling_vol(rs, ROLL_WINDOW)
rv_market <- rolling_vol(rm, ROLL_WINDOW)

# -----------------------------------------------------------------------------
# 4. CAPM REGRESSION   excess_stock = alpha + beta * excess_market + e
#    beta  = sensitivity to the market (systematic risk)
#    alpha = risk-adjusted excess return (skill / mispricing)
#    R^2   = share of the stock's moves explained by the market
# -----------------------------------------------------------------------------
excess_stock  <- rs - rf
excess_market <- rm - rf
capm <- lm(excess_stock ~ excess_market)

beta  <- coef(capm)[["excess_market"]]
alpha_weekly <- coef(capm)[["(Intercept)"]]
alpha_annual <- alpha_weekly * PERIODS_PER_YEAR
r_squared <- summary(capm)$r.squared

# -----------------------------------------------------------------------------
# 5. OUTPUT  -  summary + charts
# -----------------------------------------------------------------------------
cat("\n================ CAPM & VOLATILITY SUMMARY ================\n")
cat(sprintf("Observations (weeks)     : %d\n", length(rs)))
cat(sprintf("Beta (vs FTSE 250)       : %.2f\n", beta))
cat(sprintf("Alpha (annualised)       : %+.1f%%\n", alpha_annual * 100))
cat(sprintf("R-squared                : %.2f\n", r_squared))
cat(sprintf("Annualised vol - stock   : %.1f%%\n", vol_stock * 100))
cat(sprintf("Annualised vol - market  : %.1f%%\n", vol_market * 100))
cat(sprintf("Return correlation       : %.2f\n", cor(rs, rm)))
cat("==========================================================\n")
cat(sprintf("Interpretation: beta of %.2f means the stock is LESS volatile than\n", beta))
cat("the market (defensive). Positive alpha = outperformed its CAPM-required return.\n\n")

# Chart 1: rolling volatility
png("capm_rolling_vol.png", width = 900, height = 500, res = 130)
plot(rv_stock * 100, type = "l", lwd = 2, col = "#1f6f6b",
     main = "Rolling 26-week annualised volatility",
     xlab = "Weeks", ylab = "Annualised volatility (%)",
     ylim = range(c(rv_stock, rv_market) * 100, na.rm = TRUE))
lines(rv_market * 100, lwd = 2, col = "#2b5fb3")
legend("topright", c("A.G. Barr (BAG.L)", "FTSE 250"), col = c("#1f6f6b", "#2b5fb3"), lwd = 2, bty = "n")
dev.off()

# Chart 2: CAPM regression scatter
png("capm_scatter.png", width = 900, height = 500, res = 130)
plot(excess_market * 100, excess_stock * 100, pch = 19, col = "#2b5fb399",
     main = "CAPM regression: BAG.L excess return vs market",
     xlab = "Market excess weekly return (%)", ylab = "BAG.L excess weekly return (%)")
abline(capm, col = "#1f6f6b", lwd = 2.2)
abline(h = 0, v = 0, col = "#7a869c", lwd = 0.6)
legend("topleft", bty = "n",
       legend = c(sprintf("Beta = %.2f", beta),
                  sprintf("Alpha = %+.1f%% p.a.", alpha_annual * 100),
                  sprintf("R2 = %.2f", r_squared)))
dev.off()

cat("Charts written: capm_rolling_vol.png, capm_scatter.png\n")
