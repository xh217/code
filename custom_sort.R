# Custom sorting function
custom_sort <- function(x) {
num_part <- as.integer(gsub("\\D+", "", x))  # Extract numeric part of the strings
return(num_part)
}
