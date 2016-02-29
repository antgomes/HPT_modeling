closestinterval2 <- function( refSet, querySet, out_type=0, threshold=Inf, many_to_many=FALSE) {
  #I HAVE TO CHANGE THIS CODE TO MAKE IT MORE EFFICIENT IN R.
#I copied this one from my matlab library.
#[match] = closestinterval2(refSet,querySet,out_type,threshold,many_to_many)
#
#Same as closestinterval, but output order is switched from [querySet refSet]
#to [refSet querySet]. [refSet querySet] is more intuitive. I didn't
#change it in closestinterval because of conflict with previous scripts.
#Search for the closest interval in a reference set `RefSet` for a `querySet`.
# 
# out_type: defines the way results are output. <It would be more intuitive if
# queryset and refset output were switched>.
#    0: [refSet_index queryset_index dist]
#    1: [ref_start query_start  dis]
#    2: [ref_start ref_end query_start query_end dis]
#    3: [ref_index queryset_start dist];
#    4: [ref_start query_index dist];
# if(out_type < 0), return relative distance, otherwise, absolute distance.
# 
# threshold: 
#    if input,only queries where distance is <= threshold is returned.
# 
# many_to_many:
#    if many_to_many=true, it returns all output combination in case of
#    multiple matches.
# 
# refSet mx2 matrix, representing a set of non-overlaping
# intervals. Row i is an interval, where [m(i,1):m(i,2)] represents start and end position.
# 
# Ref(a(i,1):a(i,2)) is a query ior a single point a.
# 
# [match] = closestinterval(refSet,querySet).
# closest are the row index from refSet with closest distance to a.
# There are two possibilities to [A:B], [a,b] be overlaping intervals.
# A<a<B or a<A<b.  
# if A<a<B => A<a<b and a<B => A<b and a<B
# if a<A<b => a<A<B and A<b => A<b and a<B
# Thus:  [A:B], [a,b] are overlaping intervals <=> a<B and A<b. 
  
  if( is.null( dim(querySet) ) ){
    #If querySet is a vector, it will be turned into a matrix
    querySet =cbind(querySet, querySet);
  }
  
  if( is.null( dim(refSet) ) ){
    refSet = cbind(refSet, refSet);
  }
  
  match = data.frame(matrix(nrow=dim(querySet)[1], ncol=3));
  
  for (i in 1:dim(querySet)[1]){
    w = which( (refSet[,1] < querySet[i,2]) & (refSet[,2] > querySet[i,1]), arr.ind=TRUE ); #Check if [A:B], [a:b] overlaps. 

    if( length(w) != 0) { #if [A:B] [a:b] overlaps!
       closest = w;
#      c = w[,2];
      distance = 0;
    }
    else {
      distance = min( min( abs( cbind(refSet[,1] - querySet[i,2], refSet[,2] - querySet[i,1])  ) ) );
      w = which( abs( cbind(refSet[,1] - querySet[i,2], refSet[,2] - querySet[i,1])  ) == distance, arr.ind=TRUE);
      closest = w[,1];
      c = w[,2];
    }
    closest_u = unique(closest);

    if( (length(closest_u) >1) & !many_to_many){
      closest_u =min(closest_u);
    }
  
    for (j in 1:length(closest_u)){
      if( refSet[closest_u[j],1] > querySet[i,1] & out_type < 0) {
        distance = -distance;
      }
      match[i,] = c(closest[j], i, distance);
    }
  }

  out_type = abs(out_type);

  if(out_type == 1) {
    match = cbind( refSet[match[,1],1], querySet[match[,2], 1], match[,3]);
  }

  if(out_type == 2){
    match = cbind( refSet[match[,1],], querySet[match[,2],], match[,3]);
  }

  if(out_type == 3){
    match = cbind( match[,1], querySet[match[,2],1], match[,3]);
  }

  if(out_type == 4){
    match = cbind( refSet[match[,1],1], match[,2], match[,3]);  
  }

  ind_less_thres <- abs(match[, ncol(match)]) <= threshold
  match = match[ ind_less_thres,];

  return(match);
  
}