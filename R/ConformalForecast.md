## R Package: ConformalForecast

1. CVforecast(y, forecastfunction, h = 1, window = NULL, xreg = NULL, initial = 0, forward = TRUE, ...)
   
   - Time series cross-validation forecasting
   
   - (Similar to the tsCV function)

   - | `Input`            | (Similar to the tsCV function in the forecast pacakage)      |
     | ------------------ | ------------------------------------------------------------ |
     | `y`                | Univariate time series of length T                           |
     | `forecastfunction` | Function to return an object of class `forecast`. Its first argument must be a univariate time series, and it must have an argument `h` for the forecast horizon. If exogenous predictors are used, then it must also have `xreg` and `newxreg` arguments corresponding to the training and test periods. |
     | `h`                | Forecast horizon                                             |
     | `window`           | Length of the rolling window, if NULL, a rolling window will not be used |
     | `xreg`             | Exogeneous predictor variables passed to the forecast function if required |
     | `initial`          | Initial period of the time series where no cross-validation is performed |
     | `forward`          | If TRUE, will include h-step-ahead forecasts on forecast orgin T |
     | `...`              | Other arguments passed to `forecastfunction`                 |
     |                    |                                                              |
     
   - | `Ouput`  | an object of class "CVforecast"                              |
     | -------- | ------------------------------------------------------------ |
     | `method` | The name of the forecasting method as a character string     |
     | `mean`   | Numerical time series object containing the point forecast/s as a vector (if h=1) and a matrix otherwise |
     | `lower`  | Lower limits for prediction intervals (a list if multiple confidence levels) |
     | `upper`  | Upper limits for prediction intervals (a list if multiple confidence levels) |
     | `level`  | The confidence values associated with the prediction intervals |
     | `x`      | The original time series (either `object` itself or the time series used to create the model stored as `object`) |
     | `errors` | Forecast errors. That is x minus point forecasts             |
     |          |                                                              |

2. CPmodel(object, burnin, signedresidual = TRUE, ...)
   
   - e.g., WCP(), ACP(), PID(), MCP(),
   
   - Conformal forecasting using XX conformal prediction method
   
   - | `Input`  |                                                              |
     | -------- | ------------------------------------------------------------ |
     | `object` | An object of class "`CVforecast`". Usually the result of a call to CVforecast |
     | `burnin` | Length of burn-in period                                     |
     | `signedresidual` | If TRUE, will use signed residuals instead of abosulte residuals for conformal prediction |
     | `...`    | Other arguments for XX conformal prediction method           |
     |          |                                                              |
   
   - | `Output` | an object of class "XX"                                      |
     | -------- | ------------------------------------------------------------ |
     | `method` | The name of the conformal prediction method as a character string |
     | `mean`   | Numerical time series object containing the point forecast/s as a vector (if h=1) and a matrix otherwise |
     | `lower`  | Lower limits for prediction intervals (a list if multiple confidence levels) |
     | `upper`  | Upper limits for prediction intervals (a list if multiple confidence levels) |
     | `level`  | The confidence values associated with the prediction intervals |
     | `x`      | The original time series (either `object` itself or the time series used to create the model stored as `object`) |
     | `errors` | Forecast errors. That is x minus point forecasts             |
     | `...`    | Other outputs for XX method, like alpha, gamma, ...          |
     |          |                                                              |
   
3. Coverage(object, window = NULL)

4. Width(object, window = NULL)

5. plot()

   - coverage

   - width

   - time series plot with forecasts
