// Peter Attia
// CVsim.js
// 09-24-2017
// MIT License

var CVplot = function() {
  // Parse text field values entered by user
  var C = parseFloat(document.getElementById('conc').value);
  var D = parseFloat(document.getElementById('D').value);
  var etai = parseFloat(document.getElementById('etai').value);
  var etaf = parseFloat(document.getElementById('etaf').value);
  var v = parseFloat(document.getElementById('v').value);
  var alpha = parseFloat(document.getElementById('alpha').value);
  var k0 = parseFloat(document.getElementById('k0').value);
  var kc = parseFloat(document.getElementById('kc').value);

  // CONSTANTS
  var n = 1; // nu,ber of electrons per reaction
  var F = 96485;   // [=] C/mol, Faraday's constant
  var R = 8.3145; // [=] J/mol-K, ideal gas constant
  var T = 298.15;  // [=] K, temperature. Default = 298.15
  var f = F/(R*T); // [=] 1/V, normalized Faraday's constant at room temperature

  // SIMULATION VARIABLES
  var L  = 500;  // [=] number of iterations per t_k (pg 790). Default = 500
  var DM = 0.45; // [=] model diffusion coefficient (pg 788). Default = 0.45

  // DERIVED CONSTANTS
  var tk  = 2*(etai-etaf)/v;   // [=] s, characteristic exp. time (pg 790). In this case, total time of fwd and rev scans
  var Dt  = tk/L;              // [=] s, delta time (Eqn B.1.10, pg 790)
  var Dx  = Math.sqrt(D*Dt/DM);     // [=] cm, delta x (Eqn B.1.13, pg 791)
  var j   = Math.ceil(4.2*Math.sqrt(L))+5; // number of boxes (pg 792-793). If L~200, j=65

  // REVERSIBILITY PARAMETERS
  var ktk    = kc*tk;           // dimensionless kinetic parameter (Eqn B.3.7, pg 797)
  var km     = ktk/L;           // normalized dimensionless kinetic parameter (see bottom of pg 797)
  var Lambda = k0/Math.sqrt(D*f*v);  // dimensionless reversibility parameter (Eqn 6.4.4, pg. 236-239)

  // UPDATE REVERSIBILITY PARAMETERS
  document.getElementById("echemrev").value = Lambda.toExponential(3);
  document.getElementById("chemrev").value = km.toExponential(3);

  // WARNINGS
  var warning = '';
  if (km>0.1) {
    warning = warning.concat('k_1*t_k/L equals '+ km.toString() +', ' +
      'which exceeds the upper limit of 0.1 (see B&F, pg 797). Try lowering k_1\n');
  }
  if (C < 0) {
    warning = warning.concat('Concentration cannot be negative\n');
  }
  if (D < 0) {
    warning = warning.concat('Diffusion coefficient cannot be negative\n');
  }
  if (etai < etaf) {
    warning = warning.concat('Initial-final overpotential cannot be negative\n');
  }
  if (alpha < 0 || alpha > 1){
    warning = warning.concat('Alpha must range between 0 and 1\n')
  }
  if (k0 < 0) {
    warning = warning.concat('Electrochemical rate constant cannot be negative\n');
  }
  if (kc < 0) {
    warning = warning.concat('Chemical rate constant cannot be negative\n');
  }
  if (warning === "") {
    warning = 'No warnings';
  }
  document.getElementById("warnings").value = warning;

  // PRE-INITIALIZATION
  var k = math.range(0,L+1);    // time index vector
  var t = math.eval('Dt .* k',{Dt: Dt, k: k}); // [=] s, time vector
  var eta1 = math.eval('etai - v.*t',{etai: etai, v: v, t:t});
  var eta2 = math.eval('etaf + v.*t',{etaf: etaf, v: v, t:t});

  // overpotential scan, both directions
  var eta = [];
  var i = 0;
  var curr_eta = etai;
  while(curr_eta>etaf){
    eta[i] = curr_eta;
    i += 1;
    curr_eta = eta1.subset(math.index(i));
  }
  var k = 0;
  curr_eta = etaf;
  while(curr_eta<etai){
    eta[i] = curr_eta;
    i += 1;
    k += 1;
    curr_eta = eta2.subset(math.index(k));
  }
  eta[i] = curr_eta;

  var Enorm = math.eval('eta.*f',{eta: eta, f: f}); // dimensionless overpotential
  // kf [=] cm/s, fwd rate constant (pg 799)
  var kf = math.eval('k0*exp(  -alpha *n*Enorm)',{k0:k0, alpha:alpha, n:n, Enorm:Enorm});
  // kr [=] cm/s, rev rate constant (pg 799)
  var kb = math.eval('k0*exp((1-alpha)*n*Enorm)',{k0:k0, alpha:alpha, n:n, Enorm:Enorm});

  var O = ones(L+1,j,C); // [=] mol/cm^3, concentration of O
  var R = zeros(L+1,j);  // [=] mol/cm^3, concentration of O
  var JO = zeros1D(L+1); // [=] mol/cm^2-s, flux of O at the surface

  // START SIMULATION
  // i1 = time index. i2 = distance index
  for(var i1=0; i1<L; i1++) {
      // Update bulk concentrations of O and R
      for(var i2=1; i2<j-2; i2++) {
        //console.log(DM*(O[i1][i2+1]+O[i1][i2-1]-2*O[i1][i2]));
        O[i1+1][i2] = O[i1][i2] + DM*(O[i1][i2+1]+O[i1][i2-1]-2*O[i1][i2]);
        R[i1+1][i2] = R[i1][i2] + DM*(R[i1][i2+1]+R[i1][i2-1]-2*R[i1][i2]) - km*R[i1][i2];
      }
      // Update flux
      JO[i1+1] = ( kf[i1+1]*O[i1+1][1] - kb[i1+1]*R[i1+1][1] ) / (1 + Dx/D*(kf[i1+1] + kb[i1+1]) );
      // Update surface concentrations
      O[i1+1][0] = O[i1+1][1] - JO[i1+1]*(Dx/D);
      R[i1+1][0] = R[i1+1][1] + JO[i1+1]*(Dx/D) - km*R[i1+1][1];
  }
  // Calculate current density, Z, from flux of O
  var Z = math.eval('-n*F.*JO./10',{n:n, F:F, JO:JO}); // [=] A/m^2 -> mA/cm^2, current density

  return [eta, Z];
}

function zeros1D(rows){
  // makes a 1D array of 0s
  var arr = [];
  for (i1=0; i1<rows; i1++){
    arr[i1] = 0;
  }
  return arr;
}

function zeros(rows,cols){
  // makes a rowsxcols array of 0s
  var arr = [];
  for (i1=0; i1<rows; i1++){
    var temp = [];
    for (i2=0; i2<cols; i2++){
      temp[i2] = 0;
    }
    arr[i1] = temp;
  }
  return arr;
}


function ones(rows,cols,C){
  // not really ones. Makes a rowsxcols array of C's
  var arr = [];
  for (i1=0; i1<rows; i1++){
    var temp = [];
    for (i2=0; i2<cols; i2++){
      temp.push(C);
    }
    arr[i1] = temp;
  }
  return arr;
}
