library(survival)


wild_len <- read.csv("C:/Users/mayur/OneDrive/Documents/M2R_data/wildtype/2021-07-17WT_reporter_mother_length.csv", header=FALSE)
wild_fluor <- read.csv("C:/Users/mayur/OneDrive/Documents/M2R_data/wildtype/2021-07-17WT_reporter_mother_mean_fluor.csv", header=FALSE)
wild_time <- read.csv("C:/Users/mayur/OneDrive/Documents/M2R_data/wildtype/2021-07-17WT_reporter_time_adjusted.csv", header=FALSE)
###restrict data frame
wild_len_1 <- as.matrix(wild_len[1:189, 1:64])
wild_fluor_1 <- as.matrix(wild_fluor[1:189, 1:64])


###initialise vectors
div_len <- c()
added_len <- c()
start_1 <- c()
stop_1 <- c()
index_1 <- c()
status_1 <- c()
birth_length_1 <- c()


###update vectors
for (j in 1:64) {
  prev_stop <- 1
  prev_len_after_div <- wild_len_1[1, j]
  birth_length_1 <- c(birth_length_1, NA)
  for (i in 2:189) {
    a = wild_len_1[i-1, j]
    b = wild_len_1[i, j]
    if (!is.nan(a) && !is.nan(b)) {
      if (a-b > 5) {
        start_1 <- c(start_1, prev_stop)
        status_1 <- c(status_1, 1)
        div_len <- c(div_len, a)
        if (length(index_1) != 0) {
          if (tail(index_1, 1) == j) {
            birth_length_1 <- c(birth_length_1, prev_len_after_div)
          }
        }
        index_1 <- c(index_1, i, j)
        stop_1 <- c(stop_1, i)
        prev_stop <- i
        added_len <- c(added_len, tail(div_len, 1) - prev_len_after_div)
        prev_len_after_div <- b
      }
    }
    non_NaN_ind <- which(!is.nan(wild_len_1[, j]))    ###take care of censored
    if (i == 189){
      if (length(non_NaN_ind) != 0 && length(stop_1) != 0){
        if (tail(stop_1, 1) != tail(non_NaN_ind, 1)){
          status_1 <- c(status_1, 0)
          if (length(index_1) != 0) {
            if (tail(index_1, 1) == j) {
              birth_length_1 <- c(birth_length_1, prev_len_after_div)
            }
          }
          index_1 <- c(index_1, tail(non_NaN_ind, 1)+1, j)
          start_1 <- c(start_1, prev_stop)
          stop_1 <- c(stop_1, tail(non_NaN_ind, 1))
          div_len <- c(div_len, wild_len_1[tail(non_NaN_ind, 1), j])
          added_len <- c(added_len, tail(div_len, 1) - prev_len_after_div)
        }
      }
    }
  }
}

###want times to start at 0
start_1 <- start_1 - 1
stop_1 <- stop_1 - 1

index_1 <- matrix(index_1, ncol = 2, byrow = TRUE)
count <- as.numeric(table(index_1[, 2]))
interval_1 <- c()
for (n in count) {
  interval_1 <- c(interval_1, 1:n)
}
index_1[,1] <- index_1[,1] - 1


wt_data <- data.frame(subject = index_1[, 2], interval = interval_1, status = status_1, start = start_1, stop = stop_1, birth_length = birth_length_1, added_length = added_len, fluorescence = wild_fluor_1[index_1])
