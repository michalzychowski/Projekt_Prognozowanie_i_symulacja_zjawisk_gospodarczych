# Wczytanie bibliotek
library(dplyr)
library(zoo)
library(xts)
library(tsbox)
library(skedastic)
library(lmtest)
library(urca)
library(forecast)
library(ggplot2)

# Wczytanie danych
setwd("C:/Users/micha/Desktop/Projekt_Prognozowanie_i_symulacja_zjawisk_gospodarczych")
data <- read.csv("Data/Dane.csv", sep = ";", dec = ",", colClasses = c("character","numeric","numeric","numeric","numeric"))
data2 <- read.csv("Data/Dane2.csv", sep = ",")
data3 <- read.csv("Data/Dane3.csv", sep = ";")

## Zbiór 1
zb1 <- data.frame(data$Data,data$Podaż.pieniądza.ogółem.M3..w.mln.zł.)
colnames(zb1) <- c("Data", "Podaż pieniądza ogółem M3 w mln zł.")
head(zb1, n=10)

## Zbiór 2
zb2 <- data.frame(data$Data,data$Wskaźnik.ogólnego.klimatu.koniunktury.w.budownictwie)
colnames(zb2) <- c("Data", "Wskażnik ogólnego klimatu koniunktury w budownictwie")
head(zb2, n=10)

## Zbiór 3
zb3 <- data.frame(data$Data,data$Mieszkania.oddane.do.użytkowania..w.tys)
colnames(zb3) <- c("Data", "Mieszkania oddane do użytkowania w tys")
head(zb3, n=10)

## Zbiór 4
zb4 <- data.frame(data$Data,data$Dochody.budżetu.państwa.ogółem..w.mln.zł.)
colnames(zb4) <- c("Data", "Dochody budżetu państwa ogółem w mln zł.")
head(zb4, n=10)

## Zbiór 5
zb5 <- data.frame(data2$date,data2$meantemp)
zb5$data2.meantemp <- round(zb5$data2.meantemp, digits = 2)
colnames(zb5) <- c("Data", "Średnia temperatura w ciągu dnia.")
head(zb5, n=10)

## Zbiór 6
zb6 <- data.frame(data3$DateTime,data3$Temperature)
colnames(zb6) <- c("Data i godzina", "Temperatura.")
head(zb6, n=10)

# Szeregi czasowe
## Szereg pierwszy
ts1 <- ts(zb1$`Podaż pieniądza ogółem M3 w mln zł.`, start = c(2000, 1), frequency = 4)
plot(ts1)

## Szereg drugi
ts2 <- ts(zb2$`Wskażnik ogólnego klimatu koniunktury w budownictwie`, start = c(2000, 1), frequency = 4)
plot(ts2)

## Szereg trzeci
ts3 <- ts(zb3$`Mieszkania oddane do użytkowania w tys`, start = c(2000, 1), frequency = 4)
plot(ts3)

## Szereg czwarty
ts4 <- ts(zb4$`Dochody budżetu państwa ogółem w mln zł.`, start = c(2000, 1), frequency = 4)
plot(ts4)

## Szereg piąty
start_date <- as.Date("2017-01-01")
end_date <- as.Date("2017-04-24")
dates <- seq.Date(from = start_date, to = end_date, by = "day")
ts5 <- ts_ts(xts(zb5$`Średnia temperatura w ciągu dnia.`, order.by = dates))
plot(ts5)

## Szereg szósty
start_date2 <- as.POSIXct("2008-01-01 03:00:00")
end_date2 <- as.POSIXct("2008-12-31 23:00:00")
dates2 <- seq(from = start_date2, to = end_date2, by = "hour")
ts6 <- ts_ts(xts(zb6$Temperatura., order.by = dates2))
plot(ts6, col="blue")

# Autokorelacja
acf(ts1)
pacf(ts1)
df1 = data.frame(time=1:length(ts1),ts1)
model1 <- lm(ts1~time,data = df1)
dwtest(formula = model1, order.by = NULL)

acf(ts2)
pacf(ts2)
df2 = data.frame(time=1:length(ts2),ts2)
model2 <- lm(ts2~time,data = df2)
dwtest(formula = model2, order.by = NULL)

acf(ts3)
pacf(ts3)
df3 = data.frame(time=1:length(ts3),ts3)
model3 <- lm(ts3~time,data = df3)
dwtest(formula = model3, order.by = NULL)

acf(ts4)
pacf(ts4)
df4 = data.frame(time=1:length(ts4),ts4)
model4 <- lm(ts4~time,data = df4)
dwtest(formula = model4, order.by = NULL)

acf(ts5)
pacf(ts5)
df5 = data.frame(time=1:length(ts5),ts5)
model5 <- lm(ts5~time,data = df5)
dwtest(formula = model5, order.by = NULL)

acf(ts6)
pacf(ts6)
df6 = data.frame(time=1:length(ts6),ts6)
model6 <- lm(ts6~time,data = df6)
dwtest(formula = model6, order.by = NULL)

# Heteroskedastyczność
plot(model1$residuals)
bptest(model1)

plot(model2$residuals)
bptest(model2)

plot(model3$residuals)
bptest(model3)

white(mainlm = model3, interactions = FALSE, statonly = FALSE)
goldfeld_quandt(mainlm = model3, method = "parametric", deflator = NA, prop_central = 1/3, group1prop = 1/2, alternative = "two.sided", statonly = FALSE)

plot(model4$residuals)
bptest(model4)

plot(model5$residuals)
bptest(model5)

white(mainlm = model5, interactions = FALSE, statonly = FALSE)
goldfeld_quandt(mainlm = model5, method = "parametric", deflator = NA, prop_central = 1/3, group1prop = 1/2, alternative = "two.sided", statonly = FALSE)

plot(model6$residuals, type = "l", col="blue")
bptest(model6)

# Stacjonarność
urca::ur.kpss(ts1) %>% summary()
urca::ur.kpss(ts2) %>% summary()
urca::ur.kpss(ts3) %>% summary()
urca::ur.kpss(ts4) %>% summary()
urca::ur.kpss(ts5) %>% summary()
urca::ur.kpss(ts6) %>% summary()

# Tworzenie modeli
## Model 1
diff(ts1, lag = 4) %>% urca::ur.kpss() %>% summary()

model1$residuals %>% urca::ur.kpss() %>% summary()

ts1arima = auto.arima(ts1, seasonal = TRUE)
ts1arima

## Model 2
diff(ts2) %>% urca::ur.kpss() %>% summary()

acf(diff(ts2))
pacf(diff(ts2))

Arima(y=ts2, order = c(1,1,1), lambda = NULL)

Arima(y=ts2, seasonal = c(1,1,1), lambda = NULL)

ts2arima = Arima(y=ts2, seasonal = c(1,1,1), lambda = "auto")
ts2arima

auto.arima(ts2, seasonal = TRUE)

## Model 3

diff(ts3) %>% urca::ur.kpss() %>% summary()

acf(diff(ts3, lag = 4))
pacf(diff(ts3, lag = 4))

Arima(y=ts3, order = c(1,1,1), lambda = NULL)

Arima(y=ts3, seasonal = c(1,1,1), lambda = NULL)

ts3arima = Arima(y=ts3, seasonal = c(1,1,1), lambda = "auto")
ts3arima

auto.arima(ts3, seasonal = TRUE)

## Model 4
diff(ts4) %>% urca::ur.kpss() %>% summary()

acf(diff(ts4, lag = 4))
pacf(diff(ts4))

Arima(y=ts4, order = c(1,1,2), lambda = NULL)

Arima(y=ts4, order = c(8,2,8), lambda = NULL)

Arima(y=ts4, seasonal = c(4,2,4), lambda = NULL)

ts4arima = Arima(y=ts4, seasonal = c(4,0,4), lambda = "auto")
ts4arima

auto.arima(ts4, seasonal = TRUE)

## Model 5
diff(ts5) %>% urca::ur.kpss() %>% summary()

acf(diff(ts5))
pacf(diff(ts5))

ts5arima = Arima(y=ts5, order = c(0,1,0), lambda = NULL)
ts5arima

Arima(y=ts5, seasonal = c(0,0,0), lambda = NULL)

Arima(y=ts5, order = c(1,2,1), lambda = "auto")

auto.arima(ts5, seasonal = FALSE)

## Model 6
diff(ts6) %>% urca::ur.kpss() %>% summary()

acf(diff(ts6, lag = 24))
pacf(diff(ts6))

Arima(y=ts6, order = c(0,1,0), lambda = NULL)

Arima(y=ts6, seasonal = c(0,0,0), lambda = NULL)

Arima(y=ts6, seasonal = c(0,0,0), lambda = "auto")

ts6arima = auto.arima(ts6, seasonal = TRUE)
ts6arima

# Predykcja
## Predykcja 1
prediction1 = forecast(ts1arima, h = 12)
autoplot(ts1, series="Original Data") +
  autolayer(prediction1, series="Forecast", PI=TRUE) +
  ggtitle("Time Series and ARIMA Forecast") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_colour_manual(values=c("Original Data"="blue", "Forecast"="green"))

## Predykcja 2
prediction2 = forecast(ts2arima, h = 12)
autoplot(ts2, series="Original Data") +
  autolayer(prediction2, series="Forecast", PI=TRUE) +
  ggtitle("Time Series and ARIMA Forecast") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_colour_manual(values=c("Original Data"="blue", "Forecast"="green"))

## Predykcja 3
prediction3 = forecast(ts3arima, h = 12)
autoplot(ts3, series="Original Data") +
  autolayer(prediction3, series="Forecast", PI=TRUE) +
  ggtitle("Time Series and ARIMA Forecast") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_colour_manual(values=c("Original Data"="blue", "Forecast"="green"))

## Predykcja 4
prediction4 = forecast(ts4arima, h = 12)
autoplot(ts4, series="Original Data") +
  autolayer(prediction4, series="Forecast", PI=TRUE) +
  ggtitle("Time Series and ARIMA Forecast") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_colour_manual(values=c("Original Data"="blue", "Forecast"="green"))

## Predykcja 5
prediction5 = forecast(ts5arima, h = 12)
autoplot(ts5, series="Original Data") +
  autolayer(prediction5, series="Forecast", PI=TRUE) +
  ggtitle("Time Series and ARIMA Forecast") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_colour_manual(values=c("Original Data"="blue", "Forecast"="green"))

## Predykcja 6
prediction6 = forecast(ts6arima, h = 360)
autoplot(ts6, series="Original Data") +
  autolayer(prediction6, series="Forecast", PI=TRUE) +
  ggtitle("Time Series and ARIMA Forecast") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_colour_manual(values=c("Original Data"="blue", "Forecast"="green"))

# Dekompozycja
## Dekompozycja 1
decompose(ts1) %>% plot()

## Dekompozycja 2
decompose(ts2) %>% plot()

## Dekompozycja 3
decompose(ts3) %>% plot()

## Dekompozycja 4
decompose(ts4) %>% plot()