
calcOmega2_JW <- function(distance_matrix, metadata, ID_var, group_var) {
    # Computes the effect size omega^2 using a distance matrix for 
    #   one-way permanova test
    #
    # Args:
    #   distance_matrix: the distance matrix as class "dist";
    #   metadata: the metadata associated with the dm; metadata and dm should be 
    #       matched
    #   ID_var: the variable name that could be used to identify samples/
    #       subjects in each group; it should be quoted
    #   group_var: the variable name that stored group membership
    #
    # Returns:
    #   The calculated effect size omega^2
    
    #
    dm <- distance_matrix
    meta <- metadata
    lvls <- meta %>% pull(group_var) %>% factor() %>% levels()
    df_model <- length(lvls) - 1
    df_error <- length(as.vector(dm)) - length(lvls)
    
    ID_within <- structure(.Data = lapply(as.list(lvls), 
                                          FUN = function(x) {filter(meta, (!!as.name(group_var))==x) %>% pull(ID_var)}
    ), 
    .Names=lvls)  # ID with each level of the group variable
    dm_within <- structure(.Data = lapply(as.list(lvls), 
                                          FUN = function(x) {subset_dm_by_names(dm, ID_within[[x]])}
    ), 
    .Names=lvls)
    ss_total <- sum(as.vector(dm)^2)/attr(dm, "Size")
    ss_error <- sum(sapply(dm_within, FUN = function(x) {sum(as.vector(x)^2)/attr(x, "Size")}))
    ss_model <- ss_total - ss_error
    mss_model <- ss_model / df_model
    mss_error <- ss_error / df_error
    omega2 <- (ss_model - df_model * mss_error) / (ss_total + mss_error)
    for (i in seq(length(lvls))){
        print(paste0(lvls[i], ": ", length(ID_within[[lvls[i]]])))
    }
    return(omega2)
}