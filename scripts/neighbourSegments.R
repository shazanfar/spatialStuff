
# given the center of a polygon, expand or contract the polygon vertices 
# by a factor
expandPoint = function(x0,y0,x,y,c) {
    # x0, y0 is the centre x,y-coords
    # x, y is the point x,y-coords
    # c is the factor to multiply vector by
    
    dy = y - y0
    dx = x - x0
    tan_sign = tan(atan2(dy,dx))
    
    dc = ifelse(c >= 1, 1-c, c-1)
    csign = ifelse(c >= 1, 1, -1)
    
    if (dx == 0) {
        
        y_e = y0 + c*dy
        
        return(c(x, y_e))    
    }
    
    if (dy == 0) {
        
        x_e = x0 + c*dx
        
        return(c(x_e, y))    
    }
    
    Dx = 1 + (dy/dx)^2
    
    z2 = dx^2 + dy^2
    
    x_e = x + csign*sign(dx)*sqrt((dc^2)*z2/Dx)
    
    y_e = y + (x_e-x)*tan_sign
    
    return(c(x_e,y_e))
}

expandPointV = Vectorize(expandPoint)

closeVertices = function(x) {c(x,x[1])}



# given a dataframe - with a column to split by - output two extra columns for
# expanded vertices
addExpandedVertices = function(verticesDF,
                               xname = "segmentation_vertices_x_global_affine",
                               yname = "segmentation_vertices_y_global_affine",
                               group = "uniqueID",
                               expansionFactor = 1.1,
                               new_xname = paste0(xname, "_expanded"),
                               new_yname = paste0(yname, "_expanded")) {
    
    require(pracma)
    
    verticesDFList = split.data.frame(verticesDF,
                                      verticesDF[,group])
    
    # get expanded vertices
    verticesDFList <- lapply(verticesDFList, function(verticesDF) {
        
        closedVertices_x = closeVertices(verticesDF[,xname])
        closedVertices_y = closeVertices(verticesDF[,yname])
        
        centre = poly_center(closedVertices_x,
                             closedVertices_y)
        
        newvertices = t(expandPointV(centre[1], centre[2],
                                     verticesDF[,xname],
                                     verticesDF[,yname],
                                     c = expansionFactor))
        
        verticesDF[,new_xname] <- newvertices[,1]
        verticesDF[,new_yname] <- newvertices[,2]
        
        return(verticesDF)
    })
    
    # plot
    verticesDF_long = do.call(rbind, verticesDFList)
    
    return(verticesDF_long)
}





# given a (named) list of vertex dataframes, output a graph of overlapping
# polygons according to an expansion factor

# polygon_df_sub = subset(polygon_df, embryo == "embryo1" & pos %in% c("Pos0", "Pos1") & z == 3)
# 
# verticesDFList = split.data.frame(polygon_df_sub,
#                                   polygon_df_sub$uniqueID)

neighbourVertices = function(verticesDFList,
                             xname = "segmentation_vertices_x_global_affine",
                             yname = "segmentation_vertices_y_global_affine",
                             expansionFactor = 1.1,
                             plot = FALSE,
                             plot2 = TRUE,
                             full = FALSE,
                             verbose = FALSE) {
    
    # full means check overlaps of polygons not just vertices,
    # this might be (much?) slower
    # note this should ideally be done within a single z-stack,
    # might be useful to do over multiple but will need to think how to do it
    
    require(pracma)
    require(igraph)
    require(gtools)
    
    # get expanded vertices
    verticesDFList <- lapply(verticesDFList, function(verticesDF) {
        
        if (verbose) {
            print(verticesDF$uniqueID[1])
            }
        
        if (plot) {
            g = ggplot(verticesDF, aes(x = get(xname),
                                       y = get(yname))) + 
                geom_polygon(fill = NA, colour = "black")
            print(g)
        }
        
        closedVertices_x = closeVertices(verticesDF[,xname])
        closedVertices_y = closeVertices(verticesDF[,yname])
        
        centre = poly_center(closedVertices_x,
                             closedVertices_y)
        
        newvertices = t(expandPointV(centre[1], centre[2],
                                     verticesDF[,xname],
                                     verticesDF[,yname],
                                     c = expansionFactor))
        
        verticesDF[,"expanded_vertices_x"] <- newvertices[,1]
        verticesDF[,"expanded_vertices_y"] <- newvertices[,2]
        
        if (plot) {
            
            g + geom_polygon(aes(x = expanded_vertices_x, y = expanded_vertices_y),
                             data = verticesDF,
                             inherit.aes = FALSE,
                             fill = NA, colour = "blue")
            
        }
        
        
        return(verticesDF)
    })
    
    # plot
    verticesDF_long = do.call(rbind, verticesDFList)
    
    if (verbose) {
        print("Expanded polygons")
    }
    
    if (plot2) {
        g = ggplot(verticesDF_long, aes(x = get(xname),
                                    y = get(yname))) + 
            geom_polygon(aes(group = uniqueID), fill = NA, colour = "black") + 
            geom_polygon(aes(x = expanded_vertices_x, y = expanded_vertices_y,
                             group = uniqueID),
                         inherit.aes = FALSE, fill = NA, colour = "blue") + 
            theme_classic() + 
            NULL
        print(g)
    }
    
    # build overlap graph
    
    # do this in a loop
    verticesDFList_nms <- mixedsort(names(verticesDFList))
    neighbourSegmentsList = as.list(verticesDFList_nms)
    names(neighbourSegmentsList) <- verticesDFList_nms
    for (segment in verticesDFList_nms) {
        
        if (verbose) {
            print(segment)
        }
        
        # ask if any cell segmentation vertices are within expanded region
        # this is fast-ish
        
        inPoly = inpolygon(verticesDF_long[,xname], 
                           verticesDF_long[,yname], 
                           verticesDFList[[segment]][,"expanded_vertices_x"],
                           verticesDFList[[segment]][,"expanded_vertices_y"],
                           boundary = TRUE)
        
        neighbourSegments = unique(verticesDF_long[inPoly,"uniqueID"])
        
        if (full) {
            # check using polygon crossings
            polycrossingList = lapply(verticesDFList, function(verticesDF) {
                polycrossing = poly_crossings(t(verticesDFList[[segment]][,c("expanded_vertices_x","expanded_vertices_y")]),
                                              t(verticesDF[,c(xname,yname)]))
                if (is.null(polycrossing)) return(FALSE)
                return(TRUE)
            })
            neighbourSegmentsFull = names(verticesDFList)[unlist(polycrossingList)]
            
            neighbourSegments <- union(neighbourSegments, neighbourSegmentsFull)
        }
        
        neighbourSegmentsList[[segment]] <- sort(setdiff(neighbourSegments,segment))
        
    }
    
    neighbourSegmentsPairs = cbind(rep(names(neighbourSegmentsList),
                                       times = unlist(lapply(neighbourSegmentsList,length))),
                                   unlist(neighbourSegmentsList))
    
    if (verbose) {
        print("Got neighbour segment pairs")
    }
    
    if (plot2) {
        
        centres = do.call(rbind,lapply(verticesDFList, function(verticesDF) {
            poly_center(closeVertices(verticesDF[,xname]),
                        closeVertices(verticesDF[,yname]))
        }))
        
        centresDF = data.frame(
            segment = rownames(centres),
            segment_x = centres[,1],
            segment_y = centres[,2]
        )
        
        neighbourSegmentsDF = data.frame(
            segment1 = neighbourSegmentsPairs[,1],
            segment2 = neighbourSegmentsPairs[,2],
            segment1_x = centres[neighbourSegmentsPairs[,1],1],
            segment1_y = centres[neighbourSegmentsPairs[,1],2],
            segment2_x = centres[neighbourSegmentsPairs[,2],1],
            segment2_y = centres[neighbourSegmentsPairs[,2],2]
        )
        
        g = ggplot(verticesDF_long, aes(x = get(xname),
                                        y = get(yname))) +
            geom_polygon(aes(group = uniqueID), fill = NA, colour = "grey", size = 0.5) + 
            geom_point(aes(x = segment_x, y = segment_y), inherit.aes = FALSE,
                       data = centresDF) +
            geom_segment(aes(x = segment1_x, xend = segment2_x,
                             y = segment1_y, yend = segment2_y),
                         inherit.aes = FALSE,
                         data = neighbourSegmentsDF,
                         size = 0.5, colour = "blue") +
            theme_classic() + 
            NULL
        
        print(g)
        
    }
    
    neighbourSegmentsGraph = simplify(graph.edgelist(neighbourSegmentsPairs, 
                                                     directed = FALSE))
    
    if (verbose) {
        print("Got neighbour segment graph")
    }
    
    return(neighbourSegmentsGraph)
}


# out = neighbourVertices(verticesDFList,
#                         xname = "segmentation_vertices_x_global_affine",
#                         yname = "segmentation_vertices_y_global_affine",
#                         expansionFactor = 1.1,
#                         plot = FALSE,
#                         plot2 = TRUE,
#                         full = FALSE)
