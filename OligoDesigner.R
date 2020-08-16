#Assign a 15 nucleotide target sequence in lowercase atcg format, using spaces between nucleotides
target <- "c g g c g a c g t g a a c g g"

#convert the target to vector and check for validity
check_target <- function(target) {
  target_v <- unlist(strsplit(target, split = " "))
  print(target_v)
  if (target_v[9] != "t") {
    print("Error: the OMegazyme target must contain a U for RNA cleavage")
  } else if (target_v[8] == "c") {
    print("Warning! OMegazyme target should contain a A/G-U")
  } else if (target_v[8] == "t") {
    print("Warning! OMegazyme target should contain a A/G-U")
  } else if (length(target_v) != 15) {
    print("Error: the OMegazyme target must be 15 nucleotides in length")
  } else {
    print("The OMegazyme target is valid.")
    return(target_v)
  }
}
valid_target_v <- check_target(target)

#reverse complement the target sequence
rev_vtv <- rev(valid_target_v)
complement_a <- function(x) {
  if (x == "a") {
    sub("a", "U", x)
  } else {
    return(x)
  }}
complement_t <- function(x) {
  if (x == "t") {
    sub("t", "A", x)
  } else {
    return(x)
  }}
complement_c <- function(x) {
  if (x == "c") {
    sub("c", "G", x)
  } else {
    return(x)
  }}
complement_g <- function(x) {
  if (x == "g") {
    sub("g", "C", x)
  } else {
    return(x)
  }}
rca <- lapply(rev_vtv, complement_a)
rct <- lapply(rca, complement_t)
rcc <- lapply(rct, complement_c)
rc_target <- sapply(rcc, complement_g)
rc_target

#define left and right arm of the oligo
left_arm <- rc_target[1:7]
right_arm <- rc_target[9:15]
left_arm
right_arm

#add OMe groups to each arm
OMegaUp <- function(x) {
  paste("m", x)
}
left_arm_OMe_space <- OMegaUp(left_arm)
left_arm_OMe <- gsub(" ", "", left_arm_OMe_space)
right_arm_OMe_space <- OMegaUp(right_arm)
right_arm_OMe <- (gsub(" ", "", right_arm_OMe_space))
left <- paste(left_arm_OMe, collapse = " ")
right <- paste(right_arm_OMe, collapse = " ")
left
right

#define the catalytic core sequence
catalytic_core <- "G G C T A G C U A C A A C G A"

#combine the arms and the core
oligo <- paste(append(left, append(catalytic_core, right_arm_OMe)), collapse = " ")
#print the designed oligo
oligo