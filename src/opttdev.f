      subroutine opttdev(veg,numplt,numspc,clusid,numcls,
     +                  maxitr,minsiz,sums,numitr,
     +                  relsum,colsum,spcsum,tmpclu)
c
c* passed in/out
c
      integer numplt
      integer numspc
      integer maxitr
      integer numcls
      integer minsiz
      integer numitr
      integer clusid(numplt)
      double precision veg(numplt,numspc)
      double precision sums(0:maxitr)
      double precision relsum(numcls)
      double precision colsum(numcls)
      double precision spcsum(numspc)
c
c* local
c
      integer tmpclu(numplt)
      double precision base
      double precision totdev
      integer clstab(numcls)
      integer moveto
      integer mvefrm
      integer plttrn
      integer flag
c
c********************* one ******************************
c
      call clcdev(veg,numplt,numspc,clusid,numcls,totdev,
     +            relsum,colsum,spcsum)
      base = totdev
      numitr = maxitr
      sums(0) = totdev
      do 18 i=1,numcls
      clstab(i) = 0
   18 continue
      do 17 i=1,numplt
      clstab(clusid(i)) = clstab(clusid(i)) + 1
   17 continue
c
      do 19 k=1,maxitr
      flag = 0
        do 10 i=1,numplt
          do 11 j=1,numplt
          tmpclu(j) = clusid(j)
  11      continue
c
          do 12 j=1,numcls
          if (j .eq. clusid(i)) goto 12
          if (clstab(clusid(i)) .le. minsiz) goto 12
          tmpclu(i) = j
          call clcdev(veg,numplt,numspc,tmpclu,numcls,totdev)
          if (totdev .lt. base) then
            base = totdev
            moveto = tmpclu(i)
            mvefrm = clusid(i)
            plttrn = i
            flag = 1
          endif
   12     continue
   10   continue
      if (flag .eq. 1) then
        clusid(plttrn) = moveto
        clstab(moveto) = clstab(moveto) + 1
        clstab(mvefrm) = clstab(mvefrm) - 1
        sums(k) = base
      else
        numitr = k 
        return
      endif
   19 continue
c
      return
c
      end
