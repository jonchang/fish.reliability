#' Replicated shape landmarks for 7 fish families
#'
#' A dataset containing landmarks for 7 fish families, marked up by turkers
#' and experienced morphologists alike. Each digitizer analyzed 5 images 5
#' times each over the course of approximately an hour.
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
"fish_reliability"




#' Shape landmarks for 7 fish families
#'
#' A dataset containing landmarks for 147 species across 7 fish families.
#'
#' @format A data frame with 10,360 rows and 6 variables:
#' \describe{
#'   \item{worker}{hash identifying each digitizer}
#'   \item{mark}{the name of the landmark in question}
#'   \item{x}{the x-coordinate of this landmark}
#'   \item{y}{the y-coordinate of this landmark}
#'   \item{family}{the family that this specimen belonged to}
#'   \item{tip}{the species digitized}
#'   \item{time_taken}{the amount of time taken for a worker to submit this task}
#'   \item{worker_time_spent}{the amount of time a worker actually took to digitize the image}
#' }
"fish_families"
