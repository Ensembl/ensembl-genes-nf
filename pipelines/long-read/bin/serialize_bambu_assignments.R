assignment_cell_to_text <- function(value) {
    if (is.null(value) || length(value) == 0) return("")
    if (is.factor(value)) return(paste(as.character(value), collapse = ","))
    if (is.list(value)) {
        return(paste(vapply(value, assignment_cell_to_text, character(1)), collapse = ","))
    }
    paste(as.character(value), collapse = ",")
}

normalise_assignment_table <- function(table) {
    if (is.matrix(table)) table <- as.data.frame(table, stringsAsFactors = FALSE)
    if (!is.data.frame(table)) table <- as.data.frame(table, stringsAsFactors = FALSE)
    table[] <- lapply(table, function(column) {
        vapply(seq_along(column), function(index) assignment_cell_to_text(column[[index]]), character(1))
    })
    table
}

write_bambu_assignments <- function(assignments, output) {
    if (is.null(assignments) || length(assignments) == 0) {
        writeLines("read\ttranscript", output)
        return(invisible(NULL))
    }

    tables <- if (is.data.frame(assignments) || is.matrix(assignments)) {
        list(assignments)
    } else if (is.list(assignments)) {
        assignments
    } else {
        list(assignments)
    }
    tables <- lapply(tables, normalise_assignment_table)
    tables <- tables[vapply(tables, nrow, integer(1)) > 0]
    if (length(tables) == 0) {
        writeLines("read\ttranscript", output)
        return(invisible(NULL))
    }
    assignments <- do.call(rbind, tables)
    utils::write.table(assignments, output, sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(NULL)
}
