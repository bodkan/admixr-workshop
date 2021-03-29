library(admixr)

snp_data <- eigenstrat("subset")

result <- f4(
  W = c("French", "Dinka"), X = "Yoruba", Y = "Vindija", Z = "Chimp",
  data = snp_data
)

result
