functional_divergence <- function(X){
  # Preconditions:   X is a matrix with functional traits as columns and 
  #           samples/pixels as rows, values/traits in X should be normalized;
  #           apply normalize_traits function first. Function courtesy of 
  #           Fabian D. Schneider, code adapted from Matlab by Natasha Stavros
  
  # Postcondition: matrix of functional divergence based on 
  #           Villeger et al (2008), DOI: 10.1890/07-1206.1, 
  
  # sum by row and find indices where sum is >0
  # ind = sum( X, 2 ) > 0; ###### MATLAB CODE
  summed_rows <- rowSums(X, na.rm = TRUE)
  indices <- match(summed_rows[summed_rows > 0],summed_rows)
  
  # Find unique values of all pixels for rows where sum > 0
  # uniquePts = unique( X(ind,:), 'rows' ); ####### MATLAB CODE
  uniquePts <- uniquecombs(X[indices,],ordered=TRUE)
  
  # Find # of unique rows
  # nbUniquePts = size( uniquePts, 1 ); ##### MATLAB CODE
  number_UniquePts <- dim(uniquePts)[1]
  
  # 4 or more unique points needed
  if (number_UniquePts > 3){
    # set up all points by flattening matrix
    g <- rowMeans(uniquePts,na.rm = TRUE)
    
    # find the euclidean distance between unique points
    # NOTE: Fabian's code included g. Can't figure out how to do that here
    # dG = pdist( cat( 1, g, uniquePts ), 'euclidean' ); ### MATLAB CODE
    print("Calculating pairwise distance. This could take awhile...")
    dG <- RiemBase::rbase.pdist(riemfactory(uniquePts[!is.na(uniquePts)],name="euclidean"),parallel = TRUE)
    print("likely will throw error because you want a vector of pairwise distances and this returns a squareform")
    # pdist(uniquePts[!is.na(uniquePts)], metric = 'euclidean' ,p=2) ### ALTERNATE R IMPLEMENTATION

    dG <- dG(1:number_UniquePts )
    dG_mean <- mean(dG)
    
    deltaD <- apply(dG, MARGIN = c(1,2), FUN = function(x) 1/number_UniquePts * (x - dG_mean))
    #deltaD = sum( 1./nbUniquePts .* ( dG - dG_mean ) ) ### MATLAB CODE
    deltaAbsD <- apply(dG, MARGIN = c(1,2), FUN = function(x) 1/number_UniquePts * abs(x - dG_mean))
    # deltaAbsD = sum( 1./nbUniquePts .* abs( dG - dG_mean ) ) ### MATLAB CODE
    
    divergence <- ( deltaD + dG_mean ) / ( deltaAbsD + dG_mean )
  }
  
  return(divergence)
}

functional_evenness <- function(X){
  # Precondition: X is a matrix with functional traits as columns and samples/pixels as rows,
  #          values/traits in X should be normalized to values between 0 and 1
  #          see Fabian D Schneider et al. (2017) Nature Communications, 
  #          DOI: 10.1038/s41467-017-01530-3
  
  # Postcondition: matrix of functional evenness based on the minimum spanning tree,
  #          a calculation of Functional Evenness by Villéger, Mason, Mouillot 2008

  # sum by row and find indices where sum is >0
  #ind = sum( X, 2 ) > 0; ### MATLAB CODE
  summed_rows <- rowSums(X, na.rm = TRUE)
  indices <- match(summed_rows[summed_rows > 0],summed_rows)  
  
  # Subset data by indices
  #X = X(ind,:); ### MATLAB CODE
  #p = X; ### MATLAB CODE
  p = X[indices,]
  
  # Find number of rows
  # S = size( p, 1 ); ### MATLAB CODE
  S = dim(p)[1]
  
  # Find unique pixel values
  # pUnique = unique( p, 'rows' ); ### MATLAB CODE
  pUnique = uniquecombs(p,ordered=TRUE)
  # sUnique = size( pUnique, 1 ); ### MATLAB CODE
  sUnique = dim(pUnique)[1]
  
  # at least 3 points are needed
  # if S > 1 ### MATLAB CODE
  if (S>3){
    # if sUnique == 1 ### MATLAB CODE
    # out = 1; ### MATLAB CODE
    if (sUnique == 1){return(1)}
    else{ # else ### MATLAB CODE
      # weight matrix / distance matrix
      # temp = pdist( p, 'euclidean' ); ### MATLAB CODE
      # w = squareform(temp); ### MATLAB CODE
      print("Calculating pairwise distance. This could take awhile...")
      temp <- RiemBase::rbase.pdist(riemfactory(p,name="euclidean"),parallel = TRUE)
      
      # minimum spanning tree
      # sparseMat = sparse( w ); ### MATLAB CODE
      # [ST,pred] = graphminspantree( sparseMat, 'Method','Prim' ); ### MATLAB CODE
      
      # weighted evenness
      # EW = full( ST ); ### MATLAB CODE
      # EW = EW( EW ~= 0 ); ### MATLAB CODE
      
      # partial weighted evenness
      # PEW = EW ./ sum( EW ); ### MATLAB CODE
      
      # functional evenness index
      # sumMinPEW = 0; ### MATLAB CODE
      # for l = 1:S-1 ### MATLAB CODE
      # sumMinPEW = sumMinPEW + min( [ PEW(l) 1/(S-1) ] );### MATLAB CODE
      # end ### MATLAB CODE
      
      # out = ( sumMinPEW - 1/(S-1) ) / ( 1 - 1/(S-1) ); ### MATLAB CODE
      return(out)
    }
  }
}