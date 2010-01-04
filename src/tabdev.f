      subroutine tabdev(tab,nrow,ncol,mu,ntypes,devian,
     +                  totdev,pval,nitr,relsum,colsum,
     +                  spcsum,pclass,tclass)
c
c* passed
c
      integer nrow
      integer ncol
      integer ntypes
      double precision tab(nrow,ncol)
      integer mu(nrow)
      double precision devian(ncol)
      double precision totdev
      double precision pval(ncol)
      integer nitr
      double precision relsum(ntypes)
      double precision colsum(ntypes)
      double precision spcsum(ncol)
      integer pclass(nrow)
      integer tclass(nrow)
c
c* local
c
      double precision tmpdev
c
c* tabdev ******************** one ****************************
c
      do 10 i=1,ncol
        do 12 j=1,nrow
        colsum(mu(j)) = colsum(mu(j)) + tab(j,i)
        spcsum(i) = spcsum(i) + tab(j,i)
   12   continue
c
        do 13 j=1,ntypes
        relsum(j) = colsum(j)/spcsum(i)
        if (relsum(j) .gt. 0) then
          totdev = totdev - 2 * log(relsum(j)) * colsum(j)
          devian(i) = devian(i) - 2 * log(relsum(j)) * colsum(j)
        endif
        colsum(j) = 0.0
   13   continue
   10 continue
c
c* tabdev ***************** two ********************************
c
      do 20 i=1,ncol
      pval(i) = 0.0
        do 21 j=1,nitr
        tmpdev = 0.0
        call permute(mu,pclass,nrow,tclass)
          do 22 k=1,ntypes
          colsum(k) = 0.0
   22     continue
c
          do 23 k=1,nrow
          colsum(pclass(k)) = colsum(pclass(k)) + tab(k,i)
   23     continue
c
          do 24 k=1,ntypes
          relsum(k) = colsum(k)/spcsum(i)
          if (relsum(k) .gt. 0) then
            tmpdev = tmpdev - 2 * log(relsum(k)) * colsum(k)
          endif
   24     continue
        if (tmpdev .le. devian(i)) then
          pval(i) = pval(i) + 1
        endif
   21   continue
      pval(i) = (pval(i)+1) / (nitr+1) 
   20 continue
c
      return
c
      end