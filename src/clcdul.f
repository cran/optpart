      subroutine clcdul(veg,numplt,numspc,clusid,numcls,
     +           relfrq,relabu,indval,indcls,
     +           clstab,maxcls,sumind)
c
c* passed in
c
      integer numplt
      integer numspc
      integer numcls
      integer clusid(numplt)
      double precision veg(numplt,numspc)
      double precision relfrq(numspc,numcls)
      double precision relabu(numspc,numcls)
      double precision indval(numspc,numcls)
      double precision indcls(numspc)
      integer clstab(numcls)
      integer maxcls(numspc)
c
c* passed back
c
      double precision sumind
c
c* local 
c* 
      double precision maxval
      double precision totveg
      double precision sumrab
c
c*********************************** one *********************************
c
      do 10 i=1,numcls
      clstab(i) = 0
   10 continue
c
      do 11 i=1,numspc
        do 12 j=1,numcls 
        relabu(i,j) = 0.0
        relfrq(i,j) = 0.0  
        indval(i,j) = 0.0
   12   continue
   11 continue
c
      do 13 i=1,numplt
      clstab(clusid(i)) = clstab(clusid(i)) + 1
   13 continue
c
      sumind = 0.0
c
c****************************** two **********************************
c
      do 20 i=1,numspc
      totveg = 0
        do 21 j=1,numplt
        if (veg(j,i) .gt. 0) then
          totveg = totveg + veg(j,i)
          relabu(i,clusid(j)) = relabu(i,clusid(j)) + veg(j,i)
          relfrq(i,clusid(j)) = relfrq(i,clusid(j)) + 1
        endif
   21   continue
c
        sumrab = 0.0
        do 22 j=1,numcls
        relabu(i,j) = relabu(i,j) / clstab(j)
        sumrab = sumrab + relabu(i,j)
        relfrq(i,j) = relfrq(i,j) / clstab(j)
   22   continue
c
        maxval = 0
        do 23 j=1,numcls
        relabu(i,j) = relabu(i,j) / sumrab
        indval(i,j) = relabu(i,j) * relfrq(i,j)
        if (indval(i,j) .gt. maxval) then
          maxval = indval(i,j)
          maxcls(i) = j
        endif
   23   continue
      indcls(i) = maxval
   20 continue
c
      do 24 i=1,numspc
c       sumind = sumind + indcls(i) * clstab(maxcls(i))
        sumind = sumind + indcls(i)
   24 continue
c
      return
c
      end
