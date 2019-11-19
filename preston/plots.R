rm(list=ls())
library(tidyverse)
library(gapminder)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# extract world data
world <- ne_countries(scale = "medium", returnclass = "sf")

# create a point graph of log gdp and life exp
gapminder %>% 
    ggplot(aes(gdpPercap, lifeExp, color=year)) +
    geom_point() +
    theme_classic() +
    scale_color_distiller(palette = "Spectral") +
    scale_x_log10() +
    labs(x="GDP Per Capita", y="Life Expectancy", color="Year")

# only look at 2007 and draw a best fit line
gapminder %>% 
    filter(year == 2007) %>% 
    ggplot(aes(gdpPercap, lifeExp)) +
    geom_point(aes(color=continent)) +
    theme_classic() +
    scale_x_log10() +
    labs(x="GDP Per Capita", y="Life Expectancy", color="Continent") +
    geom_smooth(method="lm", se=F, color = "black", linetype = 3)

# extract the slope from a series of independent linear regressions over time
slopeDF <- t(sapply(sort(unique(gapminder$year)), function(y){
    # run regression from a specific year
    lm(lifeExp ~ log(gdpPercap), data = filter(gapminder, year == y)) %>%
        # run summary and extract the coefficients
        summary() %>%
        {c(.$coefficients[2,], rsq=.$r.squared)}})) %>%
    as_tibble() %>%
    # calculate upper and lower bounds
    mutate(lwr = Estimate - 1.98 * `Std. Error`) %>%
    mutate(upr = Estimate + 1.98 * `Std. Error`) %>%
    mutate(year = sort(unique(gapminder$year)))

# plot the slope with confidence intervals
slopeDF %>%
    ggplot(aes(year, Estimate, ymin = lwr, ymax = upr, color = year)) +
    geom_point() +
    geom_errorbar() +
    scale_color_viridis_c() +
    theme_classic() +
    labs(x = "Year", y = "Slope of Preston Curve") +
    ggtitle("Declining Effect of Wealth on Life Expectancy")

# group by year and show the difference between the variance of 
# gdp and life exp over time
gapminder %>%
    group_by(year) %>%
    summarize(
        `Life Expectancy (sd)` = sd(lifeExp),
        `Log GDP (sd)` = sd(log(gdpPercap))) %>%
    pivot_longer(-year) %>%
    ggplot(aes(x = year, y = value)) +
    geom_line() +
    facet_wrap(~name, scales = "free_y") +
    theme_classic() +
    labs(x = "Year", y = "Satndard Deviation") +
    ggtitle("Changes in Spread of Preston Curve Components")


gapminder %>%
    filter(year == 2007 | year == 1957) %>%
    arrange(country, year) %>%
    group_by(continent, country) %>%
    summarize(deltaLEX = diff(lifeExp), deltalGDP = diff(log(gdpPercap))) %>%
    {lm(deltaLEX ~ deltalGDP, data = .)} %>%
    summary()

# create a datset that is the difference between start and end of time in data
deltaDF <- gapminder %>%
    filter(year == 2007 | year == 1957) %>%
    arrange(country, year) %>%
    group_by(continent, country) %>%
    summarize(deltaLEX = diff(lifeExp), deltalGDP = diff(log(gdpPercap)))

deltaDF %>%
    ggplot(aes(x = deltalGDP, y = deltaLEX)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method="lm") +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_vline(xintercept = 0, linetype = 2, color = "red") +
    geom_label(
        aes(x = x, y = y, label = tx),
        fill = "#4b2e83", color = "white", fontface = "bold",
        data = tibble(x = -1, y = 30, tx = "Decrease GDP\nIncrease LEX")) +
    geom_label(
        aes(x = x, y = y, label = tx),
        fill = "#b7a57a", color = "white", fontface = "bold",
        data = tibble(x = -1, y = -5, tx = "Decrease GDP\nDecrease LEX")) +
    geom_label(
        aes(x = x, y = y, label = tx),
        fill = "#000000", color = "white", fontface = "bold",
        data = tibble(x = 3, y = -5, tx = "Increase GDP\nDecrease LEX")) +
    geom_label(
        aes(x = x, y = y, label = tx),
        fill = "#DCDCDC", color = "white", fontface = "bold",
        data = tibble(x = 3, y = 30, tx = "Increase GDP\nIncrease LEX")) +
    labs(y = "Change in Life Expectancy", x = "Change in Log GDP")

deltaDF %>%
    mutate(Delta = case_when(
        deltaLEX >= 0 & deltalGDP >= 0 ~ "Increase GDP\nIncrease LEX",
        deltaLEX < 0 & deltalGDP < 0 ~ "Decrease GDP\nDecrease LEX",
        deltaLEX >= 0 & deltalGDP < 0 ~ "Decrease GDP\nIncrease LEX",
        deltaLEX < 0 & deltalGDP >= 0 ~ "Increase GDP\nDecrease LEX",
        TRUE ~ NA_character_)) %>%
    mutate(geounit = as.character(country)) %>%
    mutate(geounit = ifelse(
        geounit == "United States", "United States of America", geounit)) %>%
    mutate(geounit = ifelse(
        geounit == "Korea, Dem. Rep.", "South Korea", geounit)) %>%
    mutate(geounit = ifelse(
        geounit == "Congo, Dem. Rep.", "Democratic Republic of the Congo", geounit)) %>%
    mutate(geounit = ifelse(
        geounit == "Congo, Rep.", "Republic of the Congo", geounit)) %>%
    {right_join(world, ., by="geounit")} %>%
    # filter(geounit != "Antarctica") %>%
    ggplot() +
    geom_sf(aes(fill=Delta)) +
    theme_void() +
    scale_fill_manual(values=c("#b7a57a", "#4b2e83", "#000000", "#DCDCDC"))
