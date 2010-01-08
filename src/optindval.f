      subroutine optindval(veg,numplt,numspc,clusid,numcls,
     +                  relfrq,relabu,indval,indcls,
     +                  clstab,
     +                  maxcls,sumind,
     +                  maxitr,minsiz,sums,numitr, 
     +                  tmpclu)
c
c* passed in/out
c
      integer numplt
      integer numspc
      integer clusid(numplt)
      integer numcls
      integer maxitr
      integer minsiz
      integer numitr
      double precision veg(numplt,numspc)
      double precision sums(0:maxitr)
c
c* passed through to clcdul
c
      double precision relfrq(numspc,numcls)
      double precision relabu(numspc,numcls)
      double precision indval(numspc,numcls)
      double precision indcls(numspc)
      double precision sumind
      integer clstab(numcls)
      integer maxcls(numspc)
c
c* scratch arrays
c
      integer tmpclu(numplt)

c* local
c
      double precision base
      integer moveto
      integer mvefrm
      integer plttrn
      integer flag
c
c********************* one ******************************
c
      call clcdul(veg,numplt,numspc,clusid,numcls,
     +            relfrq,relabu,indval,indcls,
     +            clstab,
     +            maxcls,sumind)
      base = sumind
      numitr = maxitr
      sums(0) = sumind
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
          call clcdul(veg,numplt,numspc,tmpclu,numcls,
     +            relfrq,relabu,indval,indcls,
     +            clstab,
     +            maxcls,sumind)
          if (sumind .gt. base) then
            base = sumind
            moveto = tmpclu(i)
            mvefrm = clusid(i)
            plttrn = i
            flag = 1
          endif
   12     continue
   10   continue
c
      if (flag .eq. 1) then
        clusid(plttrn) = moveto
        clstab(moveto) = clstab(moveto) + 1
        clstab(mvefrm) = clstab(mvefrm) - 1
        sums(k) = base
        sumind = base
      else
        numitr = k 
        return
      endif
   19 continue
c
      return
c
      end
