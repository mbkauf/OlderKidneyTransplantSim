#' Calibrated coefficients for deceased donor transplant equation
#'
#' A data set containing the coefficients for the top 100 best performing
#'   parameter sets from calibration. These coefficients are for a parametric
#'   survival regression using a Weibull distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item Candidate age at listing - 65
#'   \item Gender: male (indicator: 0-1)
#'   \item Race/ethnicity: Black (indicator: 0-1)
#'   \item Race/ethnicity: Hispanic (indicator: 0-1)
#'   \item Race/ethnicity: Other People of Color (indicator: 0-1)
#'   \item Years on dialysis at listing
#'   \item Year of listing - 2010
#'   \item Peak cPRA (continuous between 0-1)
#'   \item Blood type:AB (indicator: 0-1)
#'   \item Blood type:B (indicator: 0-1)
#'   \item Blood type:O (indicator: 0-1)
#'   \item History of Diabetes (indicator: 0-1)
#'   \item History of COPD (indicator: 0-1)
#'   \item History of Peripheral Vascular Disease (indicator: 0-1)
#'   \item History of Angina/Coronary Artery Disease (indicator: 0-1)
#'   \item OPTN Region 2 (indicator: 0-1)
#'   \item OPTN Region 3 (indicator: 0-1)
#'   \item OPTN Region 4 (indicator: 0-1)
#'   \item OPTN Region 5 (indicator: 0-1)
#'   \item OPTN Region 6 (indicator: 0-1)
#'   \item OPTN Region 7 (indicator: 0-1)
#'   \item OPTN Region 8 (indicator: 0-1)
#'   \item OPTN Region 9 (indicator: 0-1)
#'   \item OPTN Region 10 (indicator: 0-1)
#'   \item OPTN Region 11 (indicator: 0-1)
#'   \item Constant
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ddtx_coef_100
#' @format matrix with 100 rows and 26 columns
NULL

#' Calibrated coefficients for living donor transplant equation
#'
#' A data set containing the coefficients for the top 100 best performing
#'   parameter sets from calibration. These coefficients are for a parametric
#'   survival regression using a Gompertz distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item Candidate age at listing - 65
#'   \item Gender: male (indicator: 0-1)
#'   \item Race/ethnicity: Black (indicator: 0-1)
#'   \item Race/ethnicity: Hispanic (indicator: 0-1)
#'   \item Race/ethnicity: Other People of Color (indicator: 0-1)
#'   \item Years on dialysis at listing
#'   \item Year of listing - 2010
#'   \item Peak cPRA (continuous between 0-1)
#'   \item Blood type:AB (indicator: 0-1)
#'   \item Blood type:B (indicator: 0-1)
#'   \item Blood type:O (indicator: 0-1)
#'   \item History of Diabetes (indicator: 0-1)
#'   \item History of COPD (indicator: 0-1)
#'   \item History of Peripheral Vascular Disease (indicator: 0-1)
#'   \item History of Angina/Coronary Artery Disease (indicator: 0-1)
#'   \item OPTN Region 2 (indicator: 0-1)
#'   \item OPTN Region 3 (indicator: 0-1)
#'   \item OPTN Region 4 (indicator: 0-1)
#'   \item OPTN Region 5 (indicator: 0-1)
#'   \item OPTN Region 6 (indicator: 0-1)
#'   \item OPTN Region 7 (indicator: 0-1)
#'   \item OPTN Region 8 (indicator: 0-1)
#'   \item OPTN Region 9 (indicator: 0-1)
#'   \item OPTN Region 10 (indicator: 0-1)
#'   \item OPTN Region 11 (indicator: 0-1)
#'   \item Constant
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ldtx_coef_100
#' @format matrix with 100 rows and 26 columns
NULL

#' Calibrated coefficients for waitlist mortality equation
#'
#' A data set containing the coefficients for the top 100 best performing
#'   parameter sets from calibration. These coefficients are for a parametric
#'   survival regression using a Weibull distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item Candidate age at listing - 65
#'   \item Gender: male (indicator: 0-1)
#'   \item Race/ethnicity: Black (indicator: 0-1)
#'   \item Race/ethnicity: Hispanic (indicator: 0-1)
#'   \item Race/ethnicity: Other People of Color (indicator: 0-1)
#'   \item Years on dialysis at listing
#'   \item Year of listing - 2010
#'   \item Peak cPRA (continuous between 0-1)
#'   \item Blood type:AB (indicator: 0-1)
#'   \item Blood type:B (indicator: 0-1)
#'   \item Blood type:O (indicator: 0-1)
#'   \item History of Diabetes (indicator: 0-1)
#'   \item History of COPD (indicator: 0-1)
#'   \item History of Peripheral Vascular Disease (indicator: 0-1)
#'   \item History of Angina/Coronary Artery Disease (indicator: 0-1)
#'   \item OPTN Region 2 (indicator: 0-1)
#'   \item OPTN Region 3 (indicator: 0-1)
#'   \item OPTN Region 4 (indicator: 0-1)
#'   \item OPTN Region 5 (indicator: 0-1)
#'   \item OPTN Region 6 (indicator: 0-1)
#'   \item OPTN Region 7 (indicator: 0-1)
#'   \item OPTN Region 8 (indicator: 0-1)
#'   \item OPTN Region 9 (indicator: 0-1)
#'   \item OPTN Region 10 (indicator: 0-1)
#'   \item OPTN Region 11 (indicator: 0-1)
#'   \item Constant
#' }
#'
#' @docType data
#' @keywords datasets
#' @name mort_coef_100
#' @format matrix with 100 rows and 26 columns
NULL

#' Calibrated coefficients for waitlist removal equation
#'
#' A data set containing the coefficients for the top 100 best performing
#'   parameter sets from calibration. These coefficients are for a parametric
#'   survival regression using a Weibull distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item Candidate age at listing - 65
#'   \item Gender: male (indicator: 0-1)
#'   \item Race/ethnicity: Black (indicator: 0-1)
#'   \item Race/ethnicity: Hispanic (indicator: 0-1)
#'   \item Race/ethnicity: Other People of Color (indicator: 0-1)
#'   \item Years on dialysis at listing
#'   \item Year of listing - 2010
#'   \item Peak cPRA (continuous between 0-1)
#'   \item Blood type:AB (indicator: 0-1)
#'   \item Blood type:B (indicator: 0-1)
#'   \item Blood type:O (indicator: 0-1)
#'   \item History of Diabetes (indicator: 0-1)
#'   \item History of COPD (indicator: 0-1)
#'   \item History of Peripheral Vascular Disease (indicator: 0-1)
#'   \item History of Angina/Coronary Artery Disease (indicator: 0-1)
#'   \item OPTN Region 2 (indicator: 0-1)
#'   \item OPTN Region 3 (indicator: 0-1)
#'   \item OPTN Region 4 (indicator: 0-1)
#'   \item OPTN Region 5 (indicator: 0-1)
#'   \item OPTN Region 6 (indicator: 0-1)
#'   \item OPTN Region 7 (indicator: 0-1)
#'   \item OPTN Region 8 (indicator: 0-1)
#'   \item OPTN Region 9 (indicator: 0-1)
#'   \item OPTN Region 10 (indicator: 0-1)
#'   \item OPTN Region 11 (indicator: 0-1)
#'   \item Constant
#' }
#'
#' @docType data
#' @keywords datasets
#' @name remove_coef_100
#' @format matrix with 100 rows and 26 columns
NULL

#' Calibrated shape parameters for deceased donor transplant equation
#'
#' A data set containing the shape parameter for the top 100 best performing
#'   parameter sets from calibration. These shape parameters are for a
#'   parametric survival regression using a Weibull distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item rho shape parameter
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ddtx_shape_100
#' @format matrix with 100 rows and 1 column
NULL

#' Calibrated shape parameters for living donor transplant equation
#'
#' A data set containing the shape parameter for the top 100 best performing
#'   parameter sets from calibration. These shape parameters are for a
#'   parametric survival regression using a Gompertz distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item gamma shape parameter
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ldtx_shape_100
#' @format matrix with 100 rows and 1 column
NULL

#' Calibrated shape parameters for waitlist mortality equation
#'
#' A data set containing the shape parameter for the top 100 best performing
#'   parameter sets from calibration. These shape parameters are for a
#'   parametric survival regression using a Weibull distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item rho shape parameter
#' }
#'
#' @docType data
#' @keywords datasets
#' @name mort_shape_100
#' @format matrix with 100 rows and 1 column
NULL

#' Calibrated shape parameters for waitlist removal equation
#'
#' A data set containing the shape parameter for the top 100 best performing
#'   parameter sets from calibration. These shape parameters are for a
#'   parametric survival regression using a Weibull distribution.
#'   The variables for the equations are:
#'
#' \itemize{
#'   \item rho shape parameter
#' }
#'
#' @docType data
#' @keywords datasets
#' @name remove_shape_100
#' @format matrix with 100 rows and 1 column
NULL
