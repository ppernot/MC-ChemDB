# Misc functions ####
getSpecies = function(
  file = 'PhotoScheme.dat'
) {
  nbReac = 0
  reactants = products = params = type = orig = locnum = list()

  scheme  = read.fwf(file = file, widths = rep(11, 12))
  scheme  = t(apply(scheme, 1, function(x)
    gsub(" ", "", x)))
  for (i in 1:nrow(scheme)) {
    nbReac = nbReac + 1
    terms = scheme[i, 1:2]
    reactants[[nbReac]] = terms[!is.na(terms) &
                                  terms != "" & terms != "HV"]
    terms = scheme[i, 3:6]
    products[[nbReac]]  = terms[!is.na(terms) &
                                  terms != "" & terms != "HV"]
    terms = scheme[i, 7:12]
    params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
    type[[nbReac]]      = 'photo'
    locnum[[nbReac]]    = i
    orig[[nbReac]]      = file
  }
  species = levels(as.factor(unlist(reactants)))
  nbSpecies = length(species)

  # Build loss matrix
  L = matrix(0, ncol = nbSpecies, nrow = nbReac)
  for (m in 1:nbReac) {
    reac = unlist(reactants[m])
    for (n in 1:nbSpecies) {
      search = species[n]
      L[m, n] = length(which(search == reac)) # Loss
    }
  }
  colnames(L) = species

  return(
    list(
      species = species,
      L = L,
      params = params,
      products = products
    )
  )
}

getXShdf5 = function (
  species,
  source_dir = './Leiden/',
  getxs = TRUE # Set to FALSE to get only uncF
) {
  # Misc. info
  info = read.csv(
    file = paste0(source_dir, 'cross_section_properties.csv'),
    header = TRUE,
    skip = 18,
    stringsAsFactors = FALSE
  )

  # Get out if species not available in dataset
  if (!(species %in% trimws(info$species)))
    return(NULL)

  # Species OK: proceed
  rownames(info) = trimws(info$species)
  uncCode =  trimws(info[species, 'cross_sec_unc'])
  uF = switch(
    uncCode,
    'A+' = 1.2,
    'A'  = 1.3,
    'B'  = 2  ,
    'C'  = 10
  )

  xs = list()
  xs[['uncF']] = uF

  if(getxs) {
    # Cross sections
    file = paste0(source_dir, 'cross_sections/', species, '/', species, '.hdf5')
    fh5 = hdf5r::H5File$new(file, mode = 'r')

    for (key in c('photoabsorption',
                  'photodissociation',
                  'photoionisation',
                  'wavelength')) {
      xs[[key]] = fh5[[key]]$read()
    }
  }

  return(xs)
}

downSample = function(wl,xs,reso = 1) {
  # Define new grid
  lims = round(range(wl))
  wl1 = seq(min(lims), max(lims), by = reso)
  nl1 = length(wl1)

  # Compute increments on original grid
  dwl = diff(wl)
  ## Remove intervals of null width
  sel = dwl != 0
  dwl = dwl[sel]
  sel = c(sel, FALSE)
  pI = xs[sel] * dwl
  wl = wl[sel]

  # Interpolate from [ir]regular to regular grid
  p1 = rep(NA, length(wl1) - 1)
  i = 0
  for (il in wl1[1:(nl1 - 1)]) {
    i = i + 1
    # Which points within interval
    sel = wl >= il - 0.5 * reso &
      wl <  il + 0.5 * reso
    if (sum(sel) == 0) {
      # No point in interval: increase upper limit
      sel1 = which(wl > il + 0.5 * reso)[1]
      if (is.na(sel1)) {
        # No point within full range: contribution is null
        p1[i] = 0
      } else {
        # Found some points
        if (i == 1) {
          # If first point on grid assign first value
          p1[i] = pI[sel1] / dwl[sel1]
        } else {
          # If not first point, linear interpolation from previous value
          x0 = il - reso
          v0 = p1[i - 1]
          x1 = wl[sel1]
          v1 = pI[sel1] / dwl[sel1]
          p1[i] = v0 + (v1 - v0) / (x1 - x0) * reso
        }
      }
    } else {
      # At least one point in regular interval: sum contributions
      p1[i] = sum(pI[sel]) / sum(dwl[sel])
    }
  }
  # Remove last point
  wl1 = wl1[1:(length(wl1) - 1)]

  return(list(wl = wl1, xs = p1))
}
nds = function(ns,dist) {
  command=paste("echo ",ns," '",dist,"' | ./Bin/Rnested.x")
  # quotes around dist avoid shell interpretation
  tc=textConnection(system(command,intern=T))
  liste=scan(tc,quiet=TRUE)
  close(tc)
  nleaves=liste[1]
  nlast=nleaves*ns+1
  nds=matrix(liste[2:nlast],ncol=nleaves,byrow=T)
  return(nds)
}
gamDiri = function(x,ru) {   #Eq.9 in Plessis2010
  x = x/sum(x) # Ensure normalization
  return(
    4 / ru^2 * (
      sum( x*(1-x) ) /
      sum( x * sqrt(x*(1-x)) )
    ) - 1
  )
}
hierSampleOld = function(qy, ionic, mask, ru=c(0.1,0.1,0.1), nMC=500) {
  # Nested sampling when ionic and !ionic channels present

  nc = ncol(qy); nw = nrow(qy)

  # Total BR for neutral channels
  brNI = rep(0,nw)
  for (il in 1:nw)
    brNI[il] = sum(qy[il,!ionic])

  # Uncertainty parameter for Diri distribution
  # to get a relative uncertainty of ruBR,
  # averaged over all wavelengths...
  gammaNI = mean(
    apply(
      cbind(brNI,1-brNI), 1,
      function(x) gamDiri(x,ru=ru[1])
    )
  )

  if(sum(!ionic) > 1)
    gammaN = mean(
      apply(
        qy[,!ionic], 1,
        function(x) gamDiri(x,ru=ru[2])
      )
    )

  if(sum(ionic) > 1)
    gammaI = mean(
      apply(
        qy[,ionic], 1,
        function(x) gamDiri(x,ru=ru[3])
      )
    )

  qySample = array(
    data=NA,
    dim=c(nMC,nw,nc)
  )

  for (il in 1:nw) {

    if(sum(!ionic) > 1)
      brN = qy[il, !ionic] / sum(qy[il, !ionic])

    if(sum(ionic) > 1)
      brI = qy[il, ionic]  / sum(qy[il, ionic])

    # Nested Dirichlet
    stringBR = paste0(
      'Diri(',
        brNI[il],
          ifelse(
            sum(!ionic) > 1,
            paste0('*Diri(',paste0(brN,collapse=','),';',gammaN,'),'),
            ','
          ),
        1-brNI[il],
          if(sum(ionic) > 1)
            paste0('*Diri(',paste0(brI,collapse=','),';',gammaI,')'),
        ';',gammaNI,
      ')'
    )

    # Sample by Nested.x and reorder
    qySample[, il, order(ionic)] = nds(nMC, stringBR)

    # Order samples to improve wavelength-wise continuity
    for (ic in 1:nc)
      qySample[, il, ic] = sort(qySample[, il, ic])
  }

  # Apply mask and renormalize
  for (iMC in 1:nMC) {
    tmp = qySample[iMC, , ]
    tmp[mask] = 0
    tmp = tmp / rowSums(tmp)
    qySample[iMC, , ] = tmp
  }

  return(qySample)
}
hierSample  = function(qy, ionic,ru=c(0.1,0.1,0.1), nMC=500,eps=1e-4) {
  # Nested sampling when ionic and !ionic channels present
  # *** Treat only non-zero channels ***

  nc = ncol(qy); nw = nrow(qy)
  qySample = array(
    data = 0,
    dim  = c(nMC,nw,nc)
  )

  for (il in 1:nw) {

    # BR for neutral species
    brNI = sum(qy[il,!ionic])
    if(brNI < eps)
      brNI = 0
    if(brNI > 1-eps)
      brNI = 1

    if(brNI*(1-brNI) == 0) {
      # Sample neutral or ionic channels by flat Diri
      if(brNI == 1){
        # print(qy[il,!ionic])
        qySample[ , il, !ionic] = diriSample0(qy[il,!ionic], ru[2], nMC, eps)
      } else {
        # print(qy[il,ionic])
        qySample[ , il, ionic] = diriSample0(qy[il,ionic], ru[3], nMC, eps)
      }

    } else {
      # Use a hierarchical representation
      gammaNI = gamDiri(c(brNI,1-brNI),ru=ru[1])

      sel_nzN = TRUE
      if(sum(!ionic) > 1) {
        brN = qy[il, !ionic] / sum(qy[il, !ionic])
        brN[brN <= eps] = 0
        brN = brN / sum(brN)
        sel_nzN = brN != 0
        if(sum(sel_nzN) > 1)
          gammaN = gamDiri(brN[sel_nzN],ru=ru[2])
      }

      sel_nzI = TRUE
      if(sum(ionic) > 1) {
        brI = qy[il, ionic] / sum(qy[il, ionic])
        brI[brI <= eps] = 0
        brI = brI / sum(brI)
        sel_nzI = brI != 0
        if(sum(sel_nzI) > 1)
          gammaI = gamDiri(brI[sel_nzI],ru=ru[3])
      }

      # Nested Dirichlet
      stringBR = paste0(
        'Diri(',
        brNI,
        ifelse(
          sum(!ionic) > 1 & sum(sel_nzN) > 1,
          paste0('*Diri(',paste0(brN[sel_nzN],collapse=','),';',gammaN,'),'),
          ','
        ),
        1-brNI,
        if(sum(ionic) > 1 & sum(sel_nzI) > 1)
          paste0('*Diri(',paste0(brI[sel_nzI],collapse=','),';',gammaI,')'),
        ';',gammaNI,
        ')'
      )

      # Sample by Nested.x and reorder
      io = order(ionic)
      ret_nz = c(sel_nzN,sel_nzI)
      qySample[, il, io[ret_nz]] = nds(nMC, stringBR)
    }

    # Order samples to improve wavelength-wise continuity
    for (ic in 1:nc)
      qySample[, il, ic] = sort(qySample[, il, ic])
  }

  # Apply mask and renormalize
  for (iMC in 1:nMC) {
    tmp = qySample[iMC, , ]
    # tmp[mask] = 0
    tmp = tmp / rowSums(tmp)
    qySample[iMC, , ] = tmp
  }

  return(qySample)
}
diriSample0 = function(br, ru, nMC, eps) {

  qySample = matrix(0, nrow = nMC, ncol = length(br))

  # Count non-zero channels
  br = br / sum(br)
  br[br <= eps] = 0
  br = br / sum(br)
  sel_nz = br != 0

  if( sum(sel_nz) <= 1 ) {
    # 1 channel: no uncertainty
    qySample[,sel_nz] = 1

  } else {
    gamma = gamDiri(br[sel_nz],ru)

    # Dirichlet
    stringBR = paste0(
      'Diri(',
      paste0(br[sel_nz],collapse=','),
      ';',
      gamma,
      ')'
    )

    # Sample by Nested.x
    qySample[,sel_nz] = nds(nMC, stringBR)

  }

  return(qySample)
}
diriSample = function(qy, ru=0.1, nMC=500, eps=1e-4) {
  # Nested sampling when ionic and !ionic channels present

  nc = ncol(qy); nw = nrow(qy)
  qySample = array(
    data = 0,
    dim  = c(nMC,nw,nc)
  )

  for (il in 1:nw) {

    qySample[ , il, ] = diriSample0(qy[il,], ru, nMC, eps)

    # Order samples to improve wavelength-wise continuity
    for (ic in 1:nc)
      qySample[, il, ic] = sort(qySample[, il, ic])

  }

  # Apply mask and renormalize
  for (iMC in 1:nMC) {
    tmp = qySample[iMC, , ]
    # tmp[mask] = 0
    tmp = tmp / rowSums(tmp)
    qySample[iMC, , ] = tmp
  }

  return(qySample)
}

genCH4NeutralBR = function(nMC=500,reso=1,eps=5e-3) {
  # Generate BRs for the neutral channels of CH4 photolysis,
  # according to Gans2013

  # wavl grid
  wl = seq(100, 140, by = reso)

  nc = 4; nw = length(wl)
  qySample = array(
    data = 0,
    dim  = c(nMC,nw,nc)
  )

  # Get experimental BR sets
  # Order in files : CH3, CH2a=(1CH2), CH2X=(3CH2), CH
  nskip1= 36 # Lines to skip in sample 1
  b118_t = read.csv(
    file = paste0(source_dir,'rb-118.dat'),
    sep  = ' ',
    skip = nskip1,
    header = FALSE
  )

  nskip2=128 # id
  b121_t = read.csv(
    file = paste0(source_dir,'rb-121.dat'),
    sep  = ' ',
    skip = nskip2,
    header = FALSE
  )
  if(nMC > nrow(b118_t) | nMC > nrow(b118_t))
    stop('Sample size too large')

  xB2 = c(105.8, 123.6, 136.9)
  xB1 = c(100.0, 118.2, 121.6, 140.0)
  plot(wl,wl,type='n',ylim=c(0,1))
  for(iMC in 1:nMC) {

    # B2
    yB2 = c(
      rnorm(1,0.230,0.030),
      rnorm(1,0.059,0.005),
      eps
    )
    lyB2  = log(yB2)
    lmB2 = lm(lyB2~1+xB2+I(xB2^2))
    pB2  = predict(lmB2, newdata = data.frame(xB2=wl))
    pB2  = exp(pB2)
    # lines(wl,pB2,col=4)
    # points(xB2,yB2,col='orange')

    yB1 = matrix(
      nrow = length(xB1),
      ncol = 3,
      data = c(
        nds(1,'Diun(3)'),
        as.numeric(b118_t[iMC,1:3]),
        as.numeric(b121_t[iMC,1:3]),
        c(nds(1,'Diun(2)'),eps)
      ),
      byrow=TRUE
    )
    yB1 = yB1 / rowSums(yB1)
    # matpoints(xB1,yB1,col=1:3,pch=16)
    lyB1 = t(apply(yB1,1,function(x) log(x/prod(x)^(1/length(x))) ))
    pB1 = matrix(NA,nrow=nw,ncol=ncol(lyB1))
    for(i in 1:ncol(lyB1)) {
      lmB1 = lm(lyB1[,i]~1+xB1+I(xB1^2)+I(xB1^3))
      pB1[,i]  = predict(lmB1, newdata = data.frame(xB1=wl))
    }
    pB1  = t(apply(pB1,1,function(x) exp(x)/sum(exp(x))))
    # matlines(wl,pB1,col=1:3,lty=1)

    qySample[iMC,,] = cbind((1-pB2)*pB1,pB2)
    matlines(wl,qySample[iMC,,],col=1:4,lty=1)
  }

  return(qySample)
}
