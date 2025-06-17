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


birth_length <- wt_data[,6]
added_length <- wt_data[,7]
fluorescence <- wt_data[,8]
start <- wt_data[,4]
stop <- wt_data[,5]
status <- wt_data[,3]

fit <- coxph(Surv(start, stop, status) ~ birth_length + added_length + fluorescence , data=wt_data)


baseline_hazard_cumulative <- basehaz(fit, centered=FALSE)

baseline_hazard_matrix <- read.csv("C:/Users/mayur/Documents/baseline_hazard_matrix.csv")
names(baseline_hazard_matrix) <- c("time", "hazard")
baseline_hazard_matrix <- baseline_hazard_matrix[, c(2, 1)]
baseline_hazard <- baseline_hazard_matrix
baseline_hazard[1,1] <- 0


birth_mat <- matrix(nrow = length(baseline_hazard$time), ncol = 64)
added_mat <- matrix(nrow = length(baseline_hazard$time), ncol = 64)
fluor_mat <- matrix(nrow = length(baseline_hazard$time), ncol = 64)
for (j in 1:64) {
  for (t in baseline_hazard$time) {
    for (k in 1:table(wt_data$subject)[[j]]) {
      sliced <- wt_data[wt_data$subject == j, ]
      if (t %in% sliced[k, 4]:sliced[k, 5]) {
        birth_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j] <- sliced[k, 6]
      }
    }
    added_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j] <- wild_len_1[t+1, j] - birth_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j]
    fluor_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j] <- wild_fluor_1[t+1, j]
  }
}

hazard_matrix <- matrix(nrow = length(baseline_hazard$time), ncol = 64)
for (j in 1:64) {
  for (t in baseline_hazard$time){
    h0 <- baseline_hazard[which(baseline_hazard == t, arr.ind = TRUE)[1], 1]
    birth <- birth_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j]
    added <- added_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j]
    fluor <- fluor_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j]
    hazard_matrix[which(baseline_hazard == t, arr.ind = TRUE)[1], j] <- h0*exp(0.0100623*birth - 0.0098946*added + 0.0010837*fluor)
  }
}


hazard_matrix_na_to_0 <- hazard_matrix
###replace NA values with 0
hazard_matrix_na_to_0[is.na(hazard_matrix_na_to_0)] <- 0
cumulative_hazard_matrix <- apply(hazard_matrix_na_to_0, 2, cumsum)


wt_data_log <- wt_data
wt_data_log[,6] <- log(wt_data[,6])
wt_data_log[,7] <- log(wt_data[,7])
wt_data_log[,8] <- log(wt_data[,8])
colnames(wt_data_log)[6] <- "birth_length_l"
colnames(wt_data_log)[7] <- "added_length_l"
colnames(wt_data_log)[8] <- "fluorescence_l"

which(is.infinite(wt_data_log$added_length_l))
wt_data_log$added_length_l[217] <- NaN
wt_data_log$added_length_l[661] <- NaN

fit_l <- coxph(Surv(start, stop, status) ~ birth_length_l + added_length_l + fluorescence_l , data=wt_data_log)

baseline_hazard_l <- basehaz(fit_l, centered=FALSE)
baseline_hazard_matrix_l <- read.csv("C:/Users/mayur/OneDrive/Documents/baseline_hazard_matrix_l.csv")
names(baseline_hazard_matrix_l) <- c("time", "hazard")
baseline_hazard_matrix_l <- baseline_hazard_matrix_l[, c(2, 1)]
baseline_hazard_l <- baseline_hazard_matrix_l
baseline_hazard_l[1,1] <- 0

birth_mat_l <- matrix(nrow = length(baseline_hazard_l$time), ncol = 64)
added_mat_l <- matrix(nrow = length(baseline_hazard_l$time), ncol = 64)
fluor_mat_l <- matrix(nrow = length(baseline_hazard_l$time), ncol = 64)
hazard_matrix_l <- matrix(nrow = length(baseline_hazard_l$time), ncol = 64)
for (j in 1:64) {
  for (t in baseline_hazard_l$time) {
    for (k in 1:table(wt_data_log$subject)[[j]]) {
      sliced <- wt_data_log[wt_data_log$subject == j, ]
      if (t %in% sliced[k, 4]:sliced[k, 5]) {
        birth_mat_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], j] <- sliced[k, 6]
      }
    }
    added_mat_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], j] <- log(wild_len_1[t+1, j] - birth_mat[which(baseline_hazard == t, arr.ind = TRUE)[1], j])
    fluor_mat_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], j] <- log(wild_fluor_1[t+1, j])
  }
}

for (j in 1:64) {
  for (t in baseline_hazard_l$time){
    h0 <- baseline_hazard_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], 1]
    birth <- birth_mat_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], j]
    added <- added_mat_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], j]
    fluor <- fluor_mat_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], j]
    hazard_matrix_l[which(baseline_hazard_l == t, arr.ind = TRUE)[1], j] <- h0*exp(0.4724*birth - 0.3092*added + 0.3308*fluor)
  }
}


###clean up the data:
cell13_h_l <- hazard_matrix_l[, 13]
cell13_t_l <- baseline_hazard_l$time
cell13_remove <- c()
for (i in 1:length(baseline_hazard_l$time)) {
  if (is.nan(cell13_h_l[i])| is.na(cell13_h_l[i])|is.infinite(cell13_h_l[i])) {
    cell13_remove <- c(cell13_remove, i)
  }
}
cell13_h_l <- cell13_h_l[-cell13_remove]
cell13_t_l <- cell13_t_l[-cell13_remove]


cell25_h_l <- hazard_matrix_l[, 25]
cell25_t_l <- baseline_hazard_l$time
cell25_remove <- c()
for (i in 1:length(baseline_hazard_l$time)) {
  if (is.nan(cell25_h_l[i])| is.na(cell25_h_l[i])|is.infinite(cell25_h_l[i])) {
    cell25_remove <- c(cell25_remove, i)
  }
}
cell25_h_l <- cell25_h_l[-cell25_remove]
cell25_t_l <- cell25_t_l[-cell25_remove]


cell64_h_l <- hazard_matrix_l[, 64]
cell64_t_l <- baseline_hazard_l$time
cell64_remove <- c()
for (i in 1:length(baseline_hazard_l$time)) {
  if (is.nan(cell64_h_l[i])| is.na(cell64_h_l[i])|is.infinite(cell64_h_l[i])) {
    cell64_remove <- c(cell64_remove, i)
  }
}
cell64_h_l <- cell64_h_l[-cell64_remove]
cell64_t_l <- cell64_t_l[-cell64_remove]


par(mfrow = c(1, 3))

plot(cell13_t_l, cell13_h_l,
     type = 'l',
     xlab = 'Time', ylab = ' Hazard',
     xaxt = 'n', yaxt = 'n', main = 'Cell 13', cex.lab = 2, cex.main = 2)

axis(1, at = pretty(cell13_t_l, n = 10), cex.axis = 1.5)
axis(2, at = pretty(cell13_h_l, n = 10), cex.axis = 1.5)   

plot(cell25_t_l, cell25_h_l,
     type = 'l',
     xlab = 'Time', ylab = 'Hazard',
     xaxt = 'n', yaxt = 'n', main = 'Cell 25', cex.lab = 2, cex.main = 2)

axis(1, at = pretty(cell25_t_l, n = 10), cex.axis = 1.5)
axis(2, at = pretty(cell25_h_l, n = 10), cex.axis = 1.5)   

plot(cell64_t_l, cell64_h_l,
     type = 'l',
     xlab = 'Time', ylab = 'Hazard',
     xaxt = 'n', yaxt = 'n', main = 'Cell 64', cex.lab = 2, cex.main = 2)

axis(1, at = pretty(cell64_t_l, n = 10), cex.axis = 1.5)
axis(2, at = pretty(cell64_h_l, n = 10), cex.axis = 1.5)   


cell13_h_l <- cumsum(cell13_h_l)
cell25_h_l <- cumsum(cell25_h_l)
cell64_h_l <- cumsum(cell64_h_l)

par(mfrow = c(1, 3))

plot(cell13_t_l, cell13_h_l,
     type = 'l',
     xlab = 'Time', ylab = ' Cumulative Hazard',
     xaxt = 'n', yaxt = 'n', main = 'Cell 13', cex.lab = 2, cex.main = 2)

axis(1, at = pretty(cell13_t_l, n = 10), cex.axis = 1.5)
axis(2, at = pretty(cell13_h_l, n = 10), cex.axis = 1.5)   

plot(cell25_t_l, cell25_h_l,
     type = 'l',
     xlab = 'Time', ylab = 'Cumulative Hazard',
     xaxt = 'n', yaxt = 'n', main = 'Cell 25', cex.lab = 2, cex.main = 2)

axis(1, at = pretty(cell25_t_l, n = 10), cex.axis = 1.5)
axis(2, at = pretty(cell25_h_l, n = 10), cex.axis = 1.5)   

plot(cell64_t_l, cell64_h_l,
     type = 'l',
     xlab = 'Time', ylab = 'Cumulative Hazard',
     xaxt = 'n', yaxt = 'n', main = 'Cell 64', cex.lab = 2, cex.main = 2)

axis(1, at = pretty(cell64_t_l, n = 10), cex.axis = 1.5)
axis(2, at = pretty(cell64_h_l, n = 10), cex.axis = 1.5)   