#' Solving continuous time algebraic Riccati equation (CARE)
#'
#' The continuous time algebraic Riccati equation (CARE) has form:
#' A'*X + P*X - XB(R^-1) B'X+Q = 0, where X is an unknown symmetric matrix.
#' This function uses eigenvector based method to return the stable solution 
#' 
#' 
#'
#' @param A matrix A, should be square
#' @param B matrix B, should be square
#' @param Q matrix Q, should be invertable
#' @param R matrix R, should be square
#' @return matrix, the solution of the CARE, if algorithm failed, return NA
#' @export
#' @examples
#' 
#' my.X <- solvecare(diag(3),diag(3),diag(3),diag(3))
#' my.X + my.X-my.X %*% my.X
#' 
#' 
solvecare <- function(A, B, Q, R){
    
    if(nrow(A)!=ncol(A)){
        stop("A must be square")
    }
    

    if(nrow(B)!=ncol(B)){
        stop("B must be square")
    }

    if(nrow(Q)!=ncol(Q)){
        stop("Q must be square")
    }

    if(nrow(R)!=ncol(R)){
        stop("R must be square")
    }


    dimP <- nrow(A)
    

    if(nrow(B)!=dimP){
        stop("dimension mismatch for B, must be same as A")
    }

    if(nrow(Q)!=dimP){
        stop("dimension mismatch for Q, must be same as A")
    }

    if(nrow(R)!=dimP){
        stop("dimension mismatch for R, must be same as A")
    }

    te <- determinant(R)$modulus
    if(is.infinite(te)){
        stop("R in not invertable")
    }
    if(te <= -1e3){
        warning("determinant of R is quite small, check result")
    }

    res <- tryCatch(CARE_ArimotoPotter_cpp(A,B,Q,R), error = function(e){e})
    if(!("matrix" %in% class(res))){
        warning("Arimoto Potter eigen method failed with message: \n")
        warning(res$message)
        res <- matrix(NA, nrow(A),nrow(A))
    }

    return(res)

}


#' Solving discrete time algebraic Riccati equation (DARE)
#'
#' The discrete time algebraic Riccati equation (DARE) has form:
#' X=A'XA-(A'XB)(R+B'XB)^-1 (B'XA)+Q, where X is an unknown symmetric matrix.
#' This function uses eigenvector based method to return the stable solution 
#' 
#' 
#'
#' @param A matrix A, should be square
#' @param B matrix B, should be square
#' @param Q matrix Q, should be square
#' @param R matrix R, should be square
#' @return matrix, the solution of the DARE, if algorithm failed, return NA
#' @export
#' @examples
#' 
#' my.X <- solvedare(diag(3),diag(3),diag(3),diag(3))
#' 
solvedare <- function(A, B, Q, R){
    
    if(nrow(A)!=ncol(A)){
        stop("A must be square")
    }
    

    if(nrow(B)!=ncol(B)){
        stop("B must be square")
    }

    if(nrow(Q)!=ncol(Q)){
        stop("Q must be square")
    }

    if(nrow(R)!=ncol(R)){
        stop("R must be square")
    }


    dimP <- nrow(A)
    

    if(nrow(B)!=dimP){
        stop("dimension mismatch for B, must be same as A")
    }

    if(nrow(Q)!=dimP){
        stop("dimension mismatch for Q, must be same as A")
    }

    if(nrow(R)!=dimP){
        stop("dimension mismatch for R, must be same as A")
    }


    res <- tryCatch(DARE_ArimotoPotter_cpp(A,B,Q,R), error = function(e){e})
    if(!("matrix" %in% class(res))){
        warning("Arimoto Potter eigen method failed with message: \n")
        warning(res$message)
        res <- matrix(NA, nrow(A),nrow(A))
    }

    return(res)

}
