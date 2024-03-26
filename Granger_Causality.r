library(R.matlab)
library(bruceR)
library(vars)
library(openxlsx)

#### GCA
all_sig <- list()
file_list <- list.files("C:/Users/ASUS/Desktop/SCZ_tACS/ROI/FC_HC/cutting/EO_HC/", pattern = "*.mat", full.names = TRUE)
for (file_path in file_list) {
    base_name <- tools::file_path_sans_ext(basename(file_path))
    data <- readMat(file_path)
    oxyData <- data$nirsdata[[1]]
    oxyData[is.nan(oxyData)] <- 0
    # VARselect(oxyData)
    vm = VAR(oxyData, p=10)
    gca = granger_causality(vm)
    cleaned_sig <- gca$result$sig.F[-c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144)]
    new_sig <- as.integer(!nzchar(cleaned_sig))
    all_sig[[base_name]] <- new_sig
}
df <- data.frame(do.call(cbind, all_sig))
df[] <- 1 - df
write.xlsx(df, "C:/Users/ASUS/Desktop/SCZ_tACS/ROI/GrangerCausality/EO_HC.xlsx", rownames = FALSE)

####  NBS
file_list <- list.files("C:/Users/ASUS/Desktop/SCZ_tACS/ROI/FC_HC/cutting/HA_HC/", pattern = "*.mat", full.names = TRUE)
setwd("C:/Users/ASUS/Desktop/SCZ_tACS/ROI/FC_HC/cutting/HA_HC/")
for (file_path in file_list) {
    base_name <- tools::file_path_sans_ext(basename(file_path))

    data <- readMat(file_path)
    oxyData <- data$nirsdata[[1]]
    oxyData[is.nan(oxyData)] <- 0
    vm = VAR(oxyData, p=10)
    gca = granger_causality(vm)
    cleaned_sig <- gca$result$sig.F[-c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144)]
    newsig <- ifelse(cleaned_sig == "   ", 0, 1)
    empty_matrix <- matrix(0, nrow = 12, ncol = 12)
    matrix_sig <- matrix(newsig, nrow = 12, ncol = 11, byrow = TRUE)
    empty_matrix[1, 2:12] <- matrix_sig[1, 1:11]
    empty_matrix[2, 1] <- matrix_sig[2, 1]
    empty_matrix[2, 3:12] <- matrix_sig[2, 2:11]
    empty_matrix[3, 1:2] <- matrix_sig[3, 1:2]
    empty_matrix[3, 4:12] <- matrix_sig[3, 3:11]
    empty_matrix[4, 1:3] <- matrix_sig[4, 1:3]
    empty_matrix[4, 5:12] <- matrix_sig[4, 4:11]
    empty_matrix[5, 1:4] <- matrix_sig[5, 1:4]
    empty_matrix[5, 6:12] <- matrix_sig[5, 5:11]
    empty_matrix[6, 1:5] <- matrix_sig[6, 1:5]
    empty_matrix[6, 7:12] <- matrix_sig[6, 6:11]
    empty_matrix[7, 1:6] <- matrix_sig[7, 1:6]
    empty_matrix[7, 8:12] <- matrix_sig[7, 7:11]
    empty_matrix[8, 1:7] <- matrix_sig[8, 1:7]
    empty_matrix[8, 9:12] <- matrix_sig[8, 8:11]
    empty_matrix[9, 1:8] <- matrix_sig[9, 1:8]
    empty_matrix[9, 10:12] <- matrix_sig[9, 9:11]
    empty_matrix[10, 1:9] <- matrix_sig[10, 1:9]
    empty_matrix[10, 11:12] <- matrix_sig[10, 10:11]
    empty_matrix[11, 1:10] <- matrix_sig[11, 1:10]
    empty_matrix[11, 12] <- matrix_sig[11, 11]
    empty_matrix[12, 1:11] <- matrix_sig[12, 1:11]

    write.table(empty_matrix, file = paste0("HA_", base_name, '.txt'), sep = " ", row.names = FALSE, col.names = FALSE)
}
