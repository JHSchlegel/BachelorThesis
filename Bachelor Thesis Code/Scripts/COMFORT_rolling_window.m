load("StockPlrets.mat")

window_length = 1000;
n_windows = length(StockPlrets)-window_length;

COMFORT_VaR_mat = zeros(n_windows, 2);
COMFORT_CVaR_mat = zeros(n_windows, 2);

for t = 1:n_windows
    COMFORT_01 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GARCH', 'TS-CCC'}}, 0, 0.01, [0.05, 0.1, 0.2], 0.5);
    COMFORT_05 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GARCH', 'TS-CCC'}}, 0, 0.05, [0.05, 0.1, 0.2], 0.5);

    COMFORT_VaR_mat(t, 1) = -COMFORT_01{1,1}.VaReqw;
    COMFORT_CVaR_mat(t, 1) = -COMFORT_01{1,1}.CVaReqw;

    COMFORT_VaR_mat(t, 2) = -COMFORT_05{1,1}.VaReqw;
    COMFORT_CVaR_mat(t, 2) = -COMFORT_05{1,1}.CVaReqw;
    disp("Completed " + t + " of " + n_windows)
end



Normal_VaR_mat = zeros(n_windows, 2);
Normal_CVaR_mat = zeros(n_windows, 2);

for t = 1:n_windows
    Normal_01 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MN', 'nonAR-S', 'TS-GARCH', 'TS-DCC'}}, 0, 0.01, [0.05, 0.1, 0.2], 0.5);
    Normal_05 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MN', 'nonAR-S', 'TS-GARCH', 'TS-DCC'}}, 0, 0.05, [0.05, 0.1, 0.2], 0.5);

    Normal_VaR_mat(t, 1) = -Normal_01{1,1}.VaReqw;
    Normal_CVaR_mat(t, 1) = -Normal_01{1,1}.CVaReqw;

    Normal_VaR_mat(t, 2) = -Normal_05{1,1}.VaReqw;
    Normal_CVaR_mat(t, 2) = -Normal_05{1,1}.CVaReqw;
    disp("Completed " + t + " of " + n_windows)
end

COMFORT_VaR_DCC_GJR_mat = zeros(n_windows, 2);
COMFORT_CVaR_DCC_GJR_mat = zeros(n_windows, 2);

for t = 1:n_windows
    COMFORT_01 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GJR', 'TS-DCC'}}, 0, 0.01, [0.05, 0.1, 0.2], 0.5);
    COMFORT_05 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GJR', 'TS-DCC'}}, 0, 0.05, [0.05, 0.1, 0.2], 0.5);

    COMFORT_VaR_DCC_GJR_mat(t, 1) = -COMFORT_01{1,1}.VaReqw;
    COMFORT_CVaR_DCC_GJR_mat(t, 1) = -COMFORT_01{1,1}.CVaReqw;

    COMFORT_VaR_DCC_GJR_mat(t, 2) = -COMFORT_05{1,1}.VaReqw;
    COMFORT_CVaR_DCC_GJR_mat(t, 2) = -COMFORT_05{1,1}.CVaReqw;
    disp("Completed " + t + " of " + n_windows)
end

writematrix(COMFORT_VaR_DCC_GJR_mat, 'COMFORT_VaR_DCC_GJR.csv')
writematrix(COMFORT_CVaR_DCC_GJR_mat, 'COMFORT_CVaR_DCC_GJR.csv')

COMFORT_VaR_DCC_sGARCH_mat = zeros(n_windows, 2);
COMFORT_CVaR_DCC_sGARCH_mat = zeros(n_windows, 2);

for t = 1:n_windows
    COMFORT_01 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GARCH', 'TS-DCC'}}, 0, 0.01, [0.05, 0.1, 0.2], 0.5);
    COMFORT_05 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GARCH', 'TS-DCC'}}, 0, 0.05, [0.05, 0.1, 0.2], 0.5);

    COMFORT_VaR_DCC_sGARCH_mat(t, 1) = -COMFORT_01{1,1}.VaReqw;
    COMFORT_CVaR_DCC_sGARCH_mat(t, 1) = -COMFORT_01{1,1}.CVaReqw;

    COMFORT_VaR_DCC_sGARCH_mat(t, 2) = -COMFORT_05{1,1}.VaReqw;
    COMFORT_CVaR_DCC_sGARCH_mat(t, 2) = -COMFORT_05{1,1}.CVaReqw;
    disp("Completed " + t + " of " + n_windows)
end

writematrix(COMFORT_VaR_DCC_sGARCH_mat, 'COMFORT_VaR_DCC_sGARCH.csv')
writematrix(COMFORT_CVaR_DCC_sGARCH_mat, 'COMFORT_CVaR_DCC_sGARCH.csv')

COMFORT_VaR_CCC_GJR_mat = zeros(n_windows, 2);
COMFORT_CVaR_CCC_GJR_mat = zeros(n_windows, 2);

for t = 1:n_windows
    COMFORT_01 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GJR', 'TS-CCC'}}, 0, 0.01, [0.05, 0.1, 0.2], 0.5);
    COMFORT_05 = CVaR_COMFORT_singleWindow(StockPlrets(t:(window_length + t - 1), :), {{'MVG', 'nonAR-S', 'TS-GJR', 'TS-CCC'}}, 0, 0.05, [0.05, 0.1, 0.2], 0.5);

    COMFORT_VaR_CCC_GJR_mat(t, 1) = -COMFORT_01{1,1}.VaReqw;
    COMFORT_CVaR_CCC_GJR_mat(t, 1) = -COMFORT_01{1,1}.CVaReqw;

    COMFORT_VaR_CCC_GJR_mat(t, 2) = -COMFORT_05{1,1}.VaReqw;
    COMFORT_CVaR_CCC_GJR_mat(t, 2) = -COMFORT_05{1,1}.CVaReqw;
    disp("Completed " + t + " of " + n_windows)
end