      subroutine clcdev(tab,nrow,ncol,mu,ntypes,totdev,
     +           relsum,colsum,spcsum)
c
c* passed in
c
      integer nrow
      integer ncol
      integer ntypes
      double precision tab(nrow,ncol)
      integer mu(nrow)
      double precision relsum(ntypes)
      double precision colsum(ntypes)
      double precision spcsum(ncol)
c
c* passed back
c
      double precision totdev
c
c* tabdev ******************** one ****************************
c
      totdev = 0.0
c
      do 10 i=1,ncol
      spcsum(i) = 0.0
        do 11 j=1,ntypes
        colsum(j) = 0.0
   11   continue
c
        do 12 j=1,nrow
        colsum(mu(j)) = colsum(mu(j)) + tab(j,i)
        spcsum(i) = spcsum(i) + tab(j,i)
   12   continue
c
        do 13 j=1,ntypes
        relsum(j) = colsum(j)/spcsum(i)
        if (relsum(j) .gt. 0) then
          totdev = totdev - 2 * log(relsum(j)) * colsum(j)
        endif
   13   continue
   10 continue
c
      return
c
      end
