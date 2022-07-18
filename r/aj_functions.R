# Functions for computing Aalen-Johansen estimates and generating plots

library(tidyverse)
library(lubridate)
library(survival)
library(matrixStats)

# Compute Aalen-Johansen estimates
aj_fit <- function(x, nd) {
    evnames <- c("Death", "Discharge")
    statenames <- c("Hospital", evnames)

    x <- x %>% mutate(eventm = if_else(is.na(eventm), "0", as.character(eventm)))
    ncovvals <- dplyr::count(nd) %>% as.numeric()
    sftidy <- vector(ncovvals, mode = "list")
    varnames <- names(nd)

    # calculate the AJ survival estimates
    for (i in 1:ncovvals) {
        tbl <- x %>% inner_join(nd %>% filter(row_number() == i), by = varnames)
        if (dplyr::count(tbl) == 0) {
            next
        }
        event <- tbl$eventm
        event <- factor(event, levels = c(0, evnames))
        time <- tbl$timem %>% as.numeric()
        sf <- survival::survfit(Surv(time, event) ~ 1)
        sftidy[[i]] <- as.data.frame(unclass(sf)[c("time", "pstate", "lower", "upper")])
        sftidy[[i]] <- sftidy[[i]] %>% bind_cols(nd %>% filter(row_number() == i))
    }
    sftidy <- do.call("rbind", sftidy)
    ajlong <- sftidy %>%
        pivot_longer(
            cols = c(num_range("pstate.", 1:length(statenames)), num_range("lower.", 1:length(statenames)), num_range("upper.", 1:length(statenames))),
            names_to = c("summary", "state"),
            names_sep = "\\.", values_to = "prob"
        ) %>%
        pivot_wider(names_from = "summary", values_from = "prob")
    ajlong$state <- as.character(factor(ajlong$state, labels = statenames))
    names(ajlong)[names(ajlong) == "pstate"] <- "val"
    return(ajlong)
}

# Plot HFR estimates
ajhfr <- function(x, mat = FALSE, flip = FALSE) {
    covars <- x %>%
        select(-c(time, state, val, lower, upper)) %>%
        names()
    vars <- rlang::syms(covars)
    if (covars[1] == "monthyear") {
        x <- x %>%
            mutate(monthyear = factor(monthyear, levels = level_order))
    }
    x <- x %>%
        filter(state == "Death") %>%
        group_by(!!!vars, state) %>%
        arrange(time) %>%
        slice(n()) %>%
        ungroup()
    if (mat == TRUE) {
        return(x)
    } else {
        p1 <- x %>%
            ggplot() +
            aes(x = eval(as.name(covars[1])), val, ymin = lower, ymax = upper, color = state) +
            geom_errorbar(size = 1.05, position = position_dodge2(preserve = "single")) +
            labs(y = "Probability", x = "", title = "Hospitalised fatality risk", caption = "Error bars indicate 95% confidence intervals for HFR.") +
            coord_cartesian(ylim = c(0, NA)) +
            theme_bw() +
            theme(text = element_text(size = 16)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            theme(legend.position = "none")
        if (flip == TRUE) {
            p1 <- p1 +
                scale_x_discrete(breaks = level_order_thinned, limits = rev) +
                coord_flip(xlim = c(0, NA)) +
                theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
        } else if (covars[1] == "monthyear") {
            p1 <- p1 + scale_x_discrete(breaks = level_order_thinned)
        } else {
            p1 <- p1 + aes(color = eval(as.name(covars[1])))
        }
        return(p1)
    }
}

# Plot LOS estimates
ajlos <- function(x, mat = FALSE, stat = "median", flip = FALSE) {
    covars <- x %>%
        select(-c(time, state, val, lower, upper)) %>%
        names()
    vars <- rlang::syms(covars)
    if (stat == "median") {
        x <- x %>%
            filter(state %in% c("Death", "Discharge")) %>%
            group_by(!!!vars, state) %>%
            arrange(state, !!!vars) %>%
            mutate(
                lower = if_else(val == 0, 0, lower),
                upper = if_else(val == 0, 0, upper),
                weight = lead(val, default = 0) - val,
                uppweight = lead(upper, default = 0) - upper,
                lowweight = lead(lower, default = 0) - lower
            ) %>%
            summarise(
                med = matrixStats::weightedMedian(time, weight),
                upper = matrixStats::weightedMedian(time, uppweight),
                lower = matrixStats::weightedMedian(time, lowweight),
                .groups = "drop"
            )
    }
    if (stat == "mean") {
        x <- x %>%
            filter(state %in% c("Death", "Discharge")) %>%
            group_by(!!!vars, state) %>%
            arrange(state, !!!vars) %>%
            mutate(
                lower = if_else(val == 0, 0, lower),
                upper = if_else(val == 0, 0, upper),
                weight = lead(val, default = 0) - val,
                uppweight = lead(upper, default = 0) - upper,
                lowweight = lead(lower, default = 0) - lower,
                weight = if_else(weight < 0, 0, weight),
                uppweight = if_else(uppweight < 0, 0, uppweight),
                lowweight = if_else(lowweight < 0, 0, lowweight)
            ) %>%
            summarise(
                med = weighted.mean(time, weight),
                upper = weighted.mean(time, uppweight),
                lower = weighted.mean(time, lowweight),
                .groups = "drop"
            )
    }
    if (covars[1] == "monthyear") {
        x <- x %>%
            mutate(monthyear = factor(monthyear, levels = level_order))
    }
    if (mat == TRUE) {
        return(x)
    } else {
        p1 <- x %>%
            ggplot() +
            aes(x = eval(as.name(covars[1])), med, color = state, min = lower, max = upper) +
            geom_point(size = 2, position = position_dodge(width = 0.3)) +
            geom_linerange(position = position_dodge(width = 0.3)) +
            labs(y = "Days", x = "", title = "Median length of stay", caption = "Line ranges indicate 95% confidence intervals around median.") +
            coord_cartesian(ylim = c(0, NA)) +
            theme_bw() +
            theme(text = element_text(size = 16)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            theme(legend.title = element_blank())
        if (flip == TRUE) {
            p1 <- p1 +
                scale_x_discrete(breaks = level_order_thinned, limits = rev) +
                coord_flip(xlim = c(0, NA)) +
                theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
        } else if (covars[1] == "monthyear") {
            p1 <- p1 + scale_x_discrete(breaks = level_order_thinned)
        }
        return(p1)
    }
}

# Join HFR and LOS estimate plots
aj_full_fit <- function(df = hosp, nd, mat = FALSE, save_name = NA, col = FALSE, flip = FALSE, fig_height = 13) {
    ajfit_df <- aj_fit(df, nd = nd)

    if (mat == TRUE) {
        hfr <- ajhfr(ajfit_df, mat = TRUE)
        los <- ajlos(ajfit_df, mat = TRUE)
        out <- list("hfr" = hfr, "los" = los)
        return(out)
    } else if (length(nd) == 2 && col == TRUE) {
        covars <- nd %>% names()
        p1 <- ajhfr(ajfit_df, flip = flip) + facet_grid(rows = vars(eval(as.name(covars[2]))))
        p2 <- ajlos(ajfit_df, flip = flip) + facet_grid(rows = vars(eval(as.name(covars[2]))))
        (p1 + p2) + plot_annotation(tag_levels = "A")
    } else if (length(nd) == 2 && col == FALSE) {
        covars <- nd %>% names()
        p1 <- ajhfr(ajfit_df) + facet_grid(~ eval(as.name(covars[2])))
        p2 <- ajlos(ajfit_df) + facet_grid(~ eval(as.name(covars[2])))
        (p1 / p2) + plot_annotation(tag_levels = "A")
    } else {
        p1 <- ajhfr(ajfit_df)
        p2 <- ajlos(ajfit_df)
        (p1 / p2) + plot_annotation(tag_levels = "A")
    }
    if (!is.na(save_name)) {
        ggsave(here(figure_path, paste0(save_name, ".pdf")), device = cairo_pdf, width = 12, height = fig_height)
    }
}
