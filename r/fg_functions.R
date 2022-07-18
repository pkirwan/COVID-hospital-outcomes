# Functions for generating Fine-Gray estimate plots

library(tidyverse)
library(lubridate)
library(survival)

fg_plot <- function(var_group, plot_title) {
    labels %>%
        full_join(
            tidy_fg_fit %>%
                bind_cols(levels %>% select(label)),
            by = "label"
        ) %>%
        mutate(
            estimate = if_else(is.na(estimate), 1, estimate),
            label = factor(label, levels = labels$label),
            ref = if_else(estimate == 1, "ref", "not ref")
        ) %>%
        filter(group %in% var_group) %>%
        ggplot() +
        aes(y = label, x = estimate, xmin = conf.low, xmax = conf.high, shape = ref, color = group) +
        geom_point(size = 1.5) +
        geom_errorbarh(height = .3) +
        scale_x_continuous(name = "Hazard ratio", trans = "log10") +
        scale_y_discrete(limits = rev, name = "") +
        geom_vline(xintercept = 1, color = "black", linetype = "dashed", alpha = .5) +
        labs(
            subtitle = plot_title,
            caption = "Shown on log scale, hollow points indicate reference group."
        ) +
        theme_minimal() +
        scale_shape_manual(values = c(19, 1)) +
        theme(legend.position = "none") +
        facet_grid(group ~ ., space = "free_y", scales = "free_y", switch = "y") +
        theme(panel.spacing = unit(1, "lines")) +
        theme(
            strip.placement = "outside",
            strip.background = element_rect(fill = NA, colour = NA),
            panel.spacing = unit(0.3, "cm"), axis.title.y = element_blank()
        ) +
        theme(text = element_text(size = 16))
}

fg_plot_vaccine <- function(var_group, plot_title) {
    labels %>%
        full_join(
            tidy_fg_fit %>%
                bind_cols(levels %>% select(label)),
            by = "label"
        ) %>%
        mutate(
            estimate = if_else(is.na(estimate), 1, estimate),
            label = factor(label, levels = labels$label),
            ref = if_else(estimate == 1, "ref", "not ref")
        ) %>%
        filter(group %in% var_group) %>%
        ggplot() +
        aes(y = label, x = estimate, xmin = conf.low, xmax = conf.high, shape = ref, color = group) +
        geom_point(size = 1.5, color = "blue") +
        geom_errorbarh(height = .3, color = "blue") +
        scale_x_continuous(name = "Hazard ratio", trans = "log10") +
        scale_y_discrete(limits = rev, name = "") +
        geom_vline(xintercept = 1, color = "black", linetype = "dashed", alpha = .5) +
        labs(
            subtitle = plot_title,
            caption = "Shown on log scale, hollow points indicate reference group."
        ) +
        theme_minimal() +
        scale_shape_manual(values = c(19, 1)) +
        theme(legend.position = "none") +
        facet_grid(group ~ ., space = "free_y", scales = "free_y", switch = "y") +
        theme(panel.spacing = unit(1, "lines")) +
        theme(
            strip.placement = "outside",
            strip.background = element_rect(fill = NA, colour = NA),
            panel.spacing = unit(0.3, "cm"), axis.title.y = element_blank()
        ) +
        theme(text = element_text(size = 16))
}

fg_plot_shift <- function(plot_title) {
    labels %>%
        full_join(
            tidy_fg_fit_shift %>%
                bind_cols(levels %>% select(label)),
            by = c("label", "shift")
        ) %>%
        distinct() %>%
        mutate(
            estimate = if_else(is.na(estimate), 1, estimate),
            ref = if_else(estimate == 1, "ref", "not ref")
        ) %>%
        ggplot() +
        aes(y = label, x = estimate, xmin = conf.low, xmax = conf.high, color = shift) +
        geom_point(size = 1.5, position = position_dodge2(width = 0.9, reverse = TRUE)) +
        geom_errorbarh(height = .3, position = position_dodge2(width = 1, reverse = TRUE)) +
        scale_x_continuous(name = "Hazard ratio", trans = "log10") +
        scale_y_discrete(limits = rev(level_order), name = "") +
        scale_color_discrete("shift") +
        geom_vline(xintercept = 1, color = "black", linetype = "dashed", alpha = .5) +
        labs(
            subtitle = plot_title,
            caption = "Shown on log scale, reference group: Jun 2020."
        ) +
        theme_minimal() +
        theme(legend.position = "bottom", legend.title = element_blank()) +
        theme(panel.spacing = unit(1, "lines")) +
        theme(
            strip.placement = "outside",
            strip.background = element_rect(fill = NA, colour = NA),
            panel.spacing = unit(0.3, "cm"), axis.title.y = element_blank()
        ) +
        theme(text = element_text(size = 16)) +
        facet_wrap(~shift, nrow = 1)
}
