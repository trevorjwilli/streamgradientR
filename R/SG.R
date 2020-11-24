
#' See if point is on Line
#'
#' @param point pair of coordinates (x,y)
#' @param line_start pair of coordinates (x,y)
#' @param line_end pair of coordinates (x,y)
#'
#' @details
#'
#' See http://rstudio-pubs-static.s3.amazonaws.com/12524_7de6eb887f2945389c5d12869388be14.html
#' for explanation and original code

pointonline <- function(point, line_start, line_end) {
  # each of input parameters is pair of coordinates [x,y]
  if (identical(point, line_start) | identical(point, line_end))
  {
    return(TRUE)
  }  # check, if the points cooincides with start/end point
  if (point[1] > max(c(line_start[1], line_end[1])) | point[1] < min(c(line_start[1],
                                                                       line_end[1])) | point[2] > max(c(line_start[2], line_end[2])) | point[2] <
      min(c(line_start[2], line_end[2]))) {
    return(FALSE)  # if the point is out of the bounding box of the line, return false
  }
  if (line_start[2] == line_end[2]) {
    slope <- 0
  } else if (line_start[1] == line_end[1]) {
    return(T)
  } else {
    slope <- (line_start[2] - line_end[2])/(line_start[1] - line_end[1])
  }
  intercept <- -slope * line_start[1] + line_start[2]
  onLine <- round(point[2], digits = 0) == round((slope * point[1] + intercept),
                                                 digits = 0)
  return(onLine)
}


#' Calculate stationing (distance) along line
#'
#' @param points A SpatialPointsDataFrame in a planar (i.e. UTM) projection
#' @param lines A SpatialLinesDataFrame in a planar (i.e. UTM) projection
#' @param maxDist Maximum distance to be input into snapPointsToLines function
#'
#' @details
#'
#' See http://rstudio-pubs-static.s3.amazonaws.com/12524_7de6eb887f2945389c5d12869388be14.html
#' for original code and explanation
#'
#' @export

calculatestationing <- function(points, lines, maxDist = NA) {

  # snap points to lines from package maptools
  snapped <- snapPointsToLines(points, lines, maxDist, withAttrs = TRUE)


  stationing <- c()

  for (i in 1:length(snapped)) {
    crds_p <- sp::coordinates(snapped[i, ])
    line <- lines[snapped[i, ]$nearest_line_id, ]
    crds_l <- sp::coordinates(line)[[1]][[1]]
    d <- 0
    for (j in 2:nrow(crds_l)) {
      onLine <- pointonline(point = crds_p, line_start = crds_l[j - 1,
                                                                ], line_end = crds_l[j, ])
      if (onLine) {
        #print("ONLINE")
        d0 <- sqrt((crds_p[1] - crds_l[j - 1, 1])^2 + (crds_p[2] - crds_l[j -
                                                                            1, 2])^2)
        stationing <- c(stationing, d + d0)
        break
      }
      d <- d + sqrt((crds_l[j, 1] - crds_l[j - 1, 1])^2 + (crds_l[j, 2] -
                                                             crds_l[j - 1, 2])^2)
      #print(d)
    }
  }

  snapped$stationing <- stationing
  return(snapped)
}

#' Updated version of maptools snapPointsToLines
#'
#' @param points A SpatialPointsDataFrame in a planar (i.e. UTM) projection
#' @param lines A SpatialLinesDataFrame in a planar (i.e. UTM) projection
#' @param maxDist Numeric value for establishing a maximum distance to avoid snapping points that are farther apart; its default value is NA.
#' @param withAttrs Boolean value for preserving (TRUE) or getting rid (FALSE) of the original point attributes. Default: TRUE. This parameter is optional.
#'
#' @details see snapPointsToLines function from maptools.
#' Updated version from 	http://rstudio-pubs-static.s3.amazonaws.com/12524_7de6eb887f2945389c5d12869388be14.html
#'
#' @export

snapPointsToLines <- function(points, lines, maxDist = NA, withAttrs = TRUE) {
  if (methods::is(points, "SpatialPoints") && missing(withAttrs))
    withAttrs = FALSE
  if (!is.na(maxDist)) {
    w = rgeos::gWithinDistance(points, lines, dist = maxDist, byid = TRUE)
    validPoints = apply(w, 2, any)
    validLines = apply(w, 1, any)
    points = points[validPoints, ]
    lines = lines[validLines, ]
  }
  d = rgeos::gDistance(points, lines, byid = TRUE)
  nearest_line_index = apply(d, 2, which.min)
  coordsLines = sp::coordinates(lines)
  coordsPoints = sp::coordinates(points)
  mNewCoords = vapply(1:length(points), function(x) maptools::nearestPointOnLine(coordsLines[[nearest_line_index[x]]][[1]],
                                                                                 coordsPoints[x, ]), FUN.VALUE = c(0, 0))
  if (!is.na(maxDist))
    nearest_line_id = as.numeric(rownames(d)[nearest_line_index]) + 1 else nearest_line_id = nearest_line_index
  if (withAttrs)
    df = cbind(points@data, nearest_line_id) else df = data.frame(nearest_line_id, row.names = names(nearest_line_index))
  sp::SpatialPointsDataFrame(coords = t(mNewCoords), data = df, proj4string = sp::CRS(sp::proj4string(points)))
}


#' Calculate segment lengths and total line length
#'
#' @param coords A matrix of x (first column) and y (second column) coordinates
#' of vertices along a line.
#'
#' Calculates the segment length of each segment along a multipart line from the
#' coordinates of the vertices of that line using the equation
#' sqrt((x2 - x1)^2 + sqrt(y2 - y1)^2)
#'
#' @export

getdists <- function(coords) {
  totdist <- c(0)
  segs <- c(0)
  d <- 0
  for(i in 2:length(coords[,1])) {
    segdist <- sqrt((coords[i,][1] - coords[i-1,][1])^2 + (coords[i,][2] - coords[i-1,][2])^2)
    segs[i] <- segdist
    d <- segdist + d
    totdist[i] <- d
  }
  out <- data.frame(cbind(coords, segs, totdist))
  colnames(out) <- c("x", "y", "seglength", "totdist")
  out
}

#' Calculates the x and y coordinates of a point along a line
#'
#' @param points A 2x4 matrix where the first two columns specify the
#' x and y values for the endpoints of a line, the third column is the,
#' segment distance for the segment where the point is the ending segment,
#' and  the fourth column specifies the distance that vertex is along the
#' whole line (0 if it is the beginning of a line)
#' @param distp The distance along the entire line of the point to be
#' calculated
#'
#' @details
#' This calculates the distance along a line using the output from
#' the getdists() function using vector algebra. See
#' https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point
#' for a mathematical explanation
#'
#' @export

calcpoint <- function(points, distp) {
  # See https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point
  # For mathematical description
  d <- points[2,3] # Get segment length
  dt <- distp - points[1,4] # get the length between the lower point and the desired point
  t <- dt/d
  xp <- (1-t)*points[1,1] + t*points[2,1]
  yp <- (1-t)*points[1,2] + t*points[2,2]
  return(c(xp, yp, distp))
}

#' Calculate the x and y values of points a specified distance
#' upstream and downstream a point on a line
#'
#' @param points A SpatialPointsDataFrame in a planar projection
#' @param lines A SpatialLinesDataFrame in a planar projection
#' @param seglength Numeric 1L The distance upstream and downstream to calculate new points
#' @param idcol Character, optional column name for an identification column from points input
#'
#' @details
#' This function calculates the x and y values that are d units away from points a known
#' distance along a lines spatial object. NOTE: it will only work if the input spatial data are in
#' a planar projection (i.e. UTM).
#'
#' @export

getupdown <- function(points, lines, seglength, idcol = NULL) {
  stations <- calculatestationing(points, lines)
  if(!is.null(idcol)) {
    stations <- dplyr::select(stations@data, idcol, "nearest_line_id", "stationing")
  } else {
    stations <- dplyr::select(stations@data, "nearest_line_id", "stationing")
    stations <- data.frame(cbind(1:nrow(stations), stations))
  }
  colnames(stations) <- c("ID", "nearest_line_id", "stationing")
  outmat <- matrix(ncol = 3)
  ID <- c()
  for(i in 1:nrow(stations)) {
    pt <- stations[i,]
    ID[i] <- as.character(pt$ID)
    line <- lines[pt$nearest_line_id,]
    dists <- getdists(sp::coordinates(line)[[1]][[1]])
    updown <- c(pt$stationing - seglength, pt$stationing + seglength)
    #if(updown[1] < 0) { # see if the point is less then the segment distance
    #  updown[1] <- 0
    #}

    #if(updown[2] > dists[nrow(dists), ]$totdist) { # See if the segment length is off the headwaters
    #
    #}
    for(j in 1:2) {
      lessthan <- which(dists$totdist < updown[j])
      morethan <- which(dists$totdist > updown[j])

      pts <- dists[c(lessthan[length(lessthan)], morethan[1]),]
      newpts <- calcpoint(pts, updown[j])

      if(length(lessthan) == 0) {
        newpts <- as.matrix(dists[1,c(1,2,4)])
        colnames(newpts) <- NULL
      }
      if(length(morethan) == 0) {
        newpts <- as.matrix(dists[nrow(dists),c(1,2,4)])
        colnames(newpts) <- NULL
      }

      outmat <- rbind(outmat, newpts)

    }

  }

  ID <- rep(ID, each = 2)
  type <- rep(c("Downstream", "Upstream"), length(ID)/2)
  #print(ID)
  #print(type)
  #print(outmat)
  out <- data.frame(ID, type, data.frame(outmat[-1,]))
  rownames(out) <- NULL
  colnames(out) <- c("ID", "Type", "x", "y", "stationing")
  out
}



