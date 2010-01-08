      subroutine optsil(simptp,clusid,numplt,numclu,maxitr,sils,
     +           numitr,simptc,pltsil,tmpclu,nabor,sumnum,sumden)
c
c* passed in/out
c
      integer numplt
      integer numclu
      integer maxitr
      integer numitr
      integer clusid(numplt)
      double precision simptp(numplt,numplt)
      double precision sils(maxitr)
      double precision sumnum(numclu)
      integer sumden(numclu)
      double precision simptc(numplt,numclu)
      double precision pltsil(numplt)
      integer nabor(numplt)
      integer tmpclu(numplt)
c
c* local
c
      double precision base
      double precision maxdel
      double precision newbas
      double precision totsil
      integer moveto
      integer plttrn
      integer flag
c
      call ptc(simptp,numplt,numclu,clusid,simptc,sumnum,sumden)
      call silho(simptc,clusid,numplt,numclu,pltsil,totsil,nabor)
      base = totsil
c
      do 19 k=1,maxitr
      maxdel = 0.0
      flag = 0
        do 10 i=1,numplt
          do 11 j=1,numplt
          tmpclu(j) = clusid(j)
  11      continue
c
          do 12 j=1,numclu
          if (j .eq. clusid(i)) goto 12
          tmpclu(i) = j
          call ptc(simptp,numplt,numclu,tmpclu,simptc,sumnum,sumden)
          call silho(simptc,tmpclu,numplt,numclu,pltsil,totsil,nabor)
          if ((totsil - base) .gt. maxdel) then
            maxdel = totsil - base
            newbas = totsil
            from = clusid(i)
            moveto = tmpclu(i)
            plttrn = i
            flag = 1
          endif
   12     continue
   10   continue
      if (flag .eq. 1) then
        clusid(plttrn) = moveto
        base = newbas
        sils(k) = base
      else
        numitr = k - 1
        return
      endif
   19 continue
c
      numitr = maxitr
c
      return
c
      end


      subroutine silho(simptc,clusid,numplt,numclu,pltsil,totsil,nabor)
c
c* passed in
c
      integer numplt
      integer numclu
      double precision simptc(numplt,numclu)
      integer clusid(numplt)
      integer nabor(numplt)
c
c* passed back
c
      double precision pltsil(numplt)
      double precision totsil 
c
c* local
c
      double precision maxsim
c
      do 10 i=1,numplt
      maxsim = 0.0
      nabor(i) = clusid(i)
        do 11 j=1,numclu
        if (j .eq. clusid(i)) goto 11
        if (simptc(i,j) .gt. maxsim) then
          maxsim = simptc(i,j)
          nabor(i) = j
        endif
   11   continue
   10 continue
c
      totsil = 0.0
      do 12 i=1,numplt
      a = 1 - simptc(i,clusid(i))
      b = 1 - simptc(i,nabor(i))
      pltsil(i) = (b - a) / max(a,b)
      totsil = totsil + pltsil(i)
   12 continue
c
      return
c
      end

      subroutine ptc(simptp,numplt,numclu,clusid,simptc,sumnum,sumden)
c
c* passed in
c
      integer numplt
      integer numclu
      double precision simptp(numplt,numplt)
      integer clusid(numplt)
      double precision sumnum(numclu)
      integer sumden(numclu)
c
c* passed back
c
      double precision simptc(numplt,numclu)
c
      do 10 i=1,numplt
        do 11 j=1,numclu
        sumnum(j) = 0.0
        sumden(j) = 0.0
   11   continue
c
        do 12 j=1,numclu
          do 13 k=1,numplt
          if (i .eq. k) goto 13
          sumnum(clusid(k)) = sumnum(clusid(k)) + simptp(i,k) 
          sumden(clusid(k)) = sumden(clusid(k)) + 1
   13     continue
   12   continue
c
        do 14 j=1,numclu
        simptc(i,j) = sumnum(j)/max(sumden(j),1)
   14   continue
   10 continue
c
      return
c
      end
