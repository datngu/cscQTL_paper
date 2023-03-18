d = fread("/Users/datn/Downloads/HW2_data.csv")

xi = 0.01
y = d$y

sum(y/length(y+xi))


