#' Shape landmarks for 7 fish families
#'
#' A dataset containing landmarks for 7 fish families, marked up by turkers
#' and experienced morphologists alike.
#'
#' @format A data frame with 124875 rows and 8 variables:
#' \describe{
#'   \item{sha1mac}{hash identifying each digitizer}
#'   \item{sequence}{an integer describing the order this image was landmarked}
#'   \item{round}{an integer describing the nth time this image was landmarked}
#'   \item{role}{a factor describing whether the digitizer was a student, morphologist, or turker}
#'   \item{mark}{the name of the landmark in question}
#'   \item{x}{the x-coordinate of this landmark}
#'   \item{y}{the y-coordinate of this landmark}
#'   \item{duration}{how long, in seconds, this image took to landmark}
#'   \item{family}{the family that this specimen belonged to}
#' }
"fish_landmarks"
