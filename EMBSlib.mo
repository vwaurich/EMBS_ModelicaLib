within ;
package EMBSlib

  model EMBS_bodyExample
    Components.EMBS_Body eMBS_Body
      annotation (Placement(transformation(extent={{-40,0},{-20,20}})));

    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
    Modelica.Mechanics.MultiBody.Forces.WorldForce worldForce(N_to_m=10)
      annotation (Placement(transformation(extent={{-20,-58},{0,-38}},
            rotation=0)));
    Modelica.Blocks.Sources.Sine sine(
      freqHz=4,
      phase=1.7,
      amplitude=0,
      offset=10)    annotation (Placement(transformation(extent={{-80,-48},{-60,
              -28}},     rotation=0)));
    Modelica.Blocks.Sources.Constant const(k=0)
      annotation (Placement(transformation(extent={{-80,-78},{-60,-58}},
            rotation=0)));
  equation
    connect(eMBS_Body.frame_ref, world.frame_b) annotation (Line(
        points={{-40,10},{-80,10}},
        color={95,95,95},
        thickness=0.5));
    connect(sine.y,worldForce. force[3]) annotation (Line(points={{-59,-38},{
            -40,-38},{-40,-46.6667},{-22,-46.6667}}, color={0,0,127}));
    connect(const.y,worldForce. force[1]) annotation (Line(points={{-59,-68},{
            -40,-68},{-40,-49.3333},{-22,-49.3333}},  color={0,0,127}));
    connect(const.y,worldForce. force[2]) annotation (Line(points={{-59,-68},{
            -40,-68},{-40,-48},{-22,-48}},  color={0,0,127}));
    connect(worldForce.frame_b, eMBS_Body.frame_node[9]) annotation (Line(
        points={{0,-48},{22,-48},{22,-50},{42,-50},{42,10},{-20,10}},
        color={95,95,95},
        thickness=0.5));
  end EMBS_bodyExample;

  model test
    parameter String SIDfileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
    parameter EMBSlib.SID_Data sid = EMBSlib.SID_Data(SIDfileName) annotation (Evaluate=true);
    parameter Real Ke_M0[3,3] = EMBSlib.ExternalFunctions.getM0( sid, "Ke", 3, 3) annotation (Evaluate=true);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
            preserveAspectRatio=false)));
  end test;

  package Components
    model EMBS_Body
      parameter Integer numNodes = 9;
      parameter Integer numModes = 3 "the number of modes given in the SID file";
      parameter String SIDfileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
      constant Integer nr0 = 3 "dimension in space";

      Real q[nq] "modal coordinates";
      Real qd[nq] = der(q);
      Real qdd[nq] = der(qd);

      Modelica.SIunits.Position r_0[3](start={0,0,0}) "position of frame of reference";
      Modelica.SIunits.Velocity v[3] = der(r_0) "velocity of frame of reference";
      Modelica.SIunits.Acceleration a[3] = der(v) "acceleration of frame of reference";

      Modelica.SIunits.AngularVelocity omega[3] = Modelica.Mechanics.MultiBody.Frames.angularVelocity2(frame_ref.R);
      Modelica.SIunits.AngularAcceleration omega_d[3] = der(omega);
      Modelica.SIunits.AngularVelocity omega_tilde[3,3] = {{0, -omega[3],omega[2]},  {omega[3],0,-omega[1]}, {-omega[2],omega[1],0}};


      //SID File Data
      parameter EMBSlib.SID_Data sid = EMBSlib.SID_Data(SIDfileName) annotation (Evaluate=true);
      //mass
      parameter Modelica.SIunits.Mass mass = EMBSlib.ExternalFunctions.getMass(sid) annotation(Evaluate=true);
      parameter Modelica.SIunits.Mass mI[nr0,nr0] = identity(nr0)*mass;

      //mdCM
      parameter Real mdCM_M0[nr0,1] = EMBSlib.ExternalFunctions.getM0( sid, "mdCM", nr0, 1) annotation (Evaluate=true);
      parameter Real mdCM_M1[nr0,nq,1] = EMBSlib.ExternalFunctions.getM1(sid,   "mdCM", nr0, nq, 1) annotation (Evaluate=true);
      Real mdCM [nr0,1] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,1,mdCM_M0,mdCM_M1,q);
      Real cm[nr0,1] = mdCM/mass "centre of mass";
      Real mdCM_tilde[nr0,nr0] = {{0, -mdCM[3,1],mdCM[2,1]},  {mdCM[3,1],0,-mdCM[1,1]}, {-mdCM[2,1],mdCM[1,1],0}};

      //J, is always in a special structure in [6,1]
      constant Integer nJ = 6 "J is always a vektor of size 6";
      parameter Real J_M0[nJ,1]= EMBSlib.ExternalFunctions.getM0( sid, "J", nJ, 1) annotation (Evaluate=true);
      parameter Real J_M1[nJ,nq,1]= EMBSlib.ExternalFunctions.getM1( sid, "J", nJ, nq, 1) annotation (Evaluate=true);
      Real J_ [nJ,1] = EMBSlib.MatrixFunctions.getTaylorFunction(nJ,nq,1,J_M0,J_M1,q);
      Real J [nr0,nr0] = {{J_[1,1],J_[4,1],J_[5,1]},{J_[4,1],J_[2,1],J_[6,1]},{J_[5,1],J_[6,1],J_[3,1]}};

      //Ct
      parameter Real Ct_M0[nq,nr0] = EMBSlib.ExternalFunctions.getM0( sid, "Ct", nq, nr0) annotation (Evaluate=true);
      parameter Real Ct_M1[nq,nq,nr0] = EMBSlib.ExternalFunctions.getM1( sid, "Ct", nq, nq, nr0) annotation (Evaluate=true);
      Real Ct [nq,nr0] = EMBSlib.MatrixFunctions.getTaylorFunction(nq,nq,nr0,Ct_M0,Ct_M1,q);


      //Cr
      parameter Real Cr_M0[nq,nr0] = EMBSlib.ExternalFunctions.getM0( sid, "Cr", nq, nr0) annotation (Evaluate=true);
      parameter Real Cr_M1[nq,nq,nr0] = EMBSlib.ExternalFunctions.getM1( sid, "Cr", nq,  nq, nr0) annotation (Evaluate=true);
      Real Cr [nq,nr0] = EMBSlib.MatrixFunctions.getTaylorFunction(nq,nq,nr0,Cr_M0,Cr_M1,q);

      //Gr
      parameter Real Gr_ [nr0,nr0*nq] = EMBSlib.ExternalFunctions.getM0( sid, "Gr", nq, nr0*nq) annotation (Evaluate=true);
      Real Gr [nr0,nr0] = EMBSlib.MatrixFunctions.getGrMatrix(nr0,nq,nr0*nq,Gr_,qd);//! this is the sum of the product of the GR matrices with the respective qd_i

      //Ge
      parameter Real Ge_ [nq,nr0*nq] = EMBSlib.ExternalFunctions.getM0( sid, "Ge", nq, nr0*nq) annotation (Evaluate=true);
      Real Ge [nq,nr0] = EMBSlib.MatrixFunctions.getGrMatrix(nq,nq,nr0*nq,Ge_,qd);//! this is the sum of the product of the Ge matrices with the respective qd_i

      //Me
      parameter Real Me[nq,nq] = identity(nq) annotation (Evaluate=true);

      //Oe
      constant Integer nOe = 6 "its always 6";
      parameter Real Oe_M0[nq,nOe] = EMBSlib.ExternalFunctions.getM0( sid, "Oe", nq, nOe) annotation (Evaluate=true);
      parameter Real Oe_M1[nq,nq,nOe] = EMBSlib.ExternalFunctions.getM1( sid, "Oe", nq,  nq, nOe) annotation (Evaluate=true);
      Real Oe_[nq,nOe] = EMBSlib.MatrixFunctions.getTaylorFunction(3,nq,6,Oe_M0,Oe_M1,q);//wallrapp says on page 12 table 1 there is acompact form, not used here?!?!?!?

      //OMega
      Real OMega[nOe] = {omega[1]^2, omega[2]^2, omega[3]^2, omega[1]*omega[2], omega[2]*omega[3], omega[1]*omega[3]};

      //ksigma
      parameter Real ksigma[nq,1] = zeros(nq,1) annotation (Evaluate=true); //!! For simplicity....change this

      //Ke
      parameter Real Ke[nq,nq] = EMBSlib.ExternalFunctions.getM0( sid, "Ke", nq, nq) annotation (Evaluate=true);


       //De
      parameter Real De[nq,nq] =  EMBSlib.ExternalFunctions.getM0( sid, "De", nq, nq) annotation (Evaluate=true);


      // kinematic equations --> translation
      Real M_t[3] = mI*a + transpose(mdCM_tilde)* omega_d + transpose(Ct)* qdd;
      Real k_omega_t[3] = 2*omega_tilde*transpose(Ct)*qd + res2_1[:,1];
      Real res2_1[3,1] = omega_tilde*omega_tilde*cm; // need the intermediate value to fix the dimensions
      //Real f_comp[3] = M_t+k_omega_t;
      Real hd_t[3] = sum(identity(3)*nodes[i].f for i in 1:numNodes);
      //Real residual_t[3] = f_comp - hd_t - frame_ref.f;  // this should be zero

      // kinematic equations --> rotation
      Real M_r[3] = mdCM_tilde*a + J* omega_d + transpose(Cr) * qdd;
      Real k_omega_r[3] = Gr* omega_d + omega_tilde*J*omega;
      //Real t_comp[3] = M_r+k_omega_r;
      Real hd_r[3] = sum(identity(3)*nodes[i].t for i in 1:numNodes);

      //kinematic equations --> modal
      Real M_q[3] = Ct*a + Cr* omega_d + Me* qdd;
      Real k_omega_q[3] =  Ge*omega + Oe_*OMega;
      Real k_q[nq] = ksigma[:,1] + Ke*q +De*qd;
      //Real q_comp[nq] = M_q + k_omega_q + k_q;
      Real hd_e[nq] = sum(nodes[i].hde_i for i in 1:numNodes);

      Node nodes[numNodes](each sid=sid, each nq=nq, nodeArrayIdx = 1:numNodes) annotation (Placement(transformation(extent={{-10,-4},
                {10,16}})));

      Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_ref
     annotation (Placement(transformation(extent={{-116,-16},{-84,16}})));
      Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_node[numNodes] annotation (Placement(transformation(extent={{84,-16},
                {116,16}})));

      Modelica.Blocks.Sources.RealExpression qExp[nq](y=q)
        annotation (Placement(transformation(extent={{-58,40},{-38,60}})));
    protected
      parameter Integer nq = numModes;

    equation

     //q = {0,0,0};
     //kinematic equations
     //M_t+k_omega_t = hd_t;
     //M_r+k_omega_r = hd_r;
     M_q + k_omega_q + k_q = hd_e;


     r_0 = frame_ref.r_0;
     //frame_ref.f = {0,0,0};
     //frame_ref.t = {0,0,0};

     for i in 1:numNodes loop
      connect(qExp.y, nodes[i].q)    annotation (Line(points={{-37,50},{-18,50},{-18,5.8},{-10.2,5.8}},  color={0,0,127}));
      connect(nodes[i].frame_a, frame_ref) annotation (Line(
          points={{-10,0.4},{-50,0.4},{-50,0},{-100,0}},
          color={95,95,95},
          thickness=0.5));
      connect(nodes[i].frame_b, frame_node[i]) annotation (Line(
          points={{10,0.2},{59,0.2},{59,0},{100,0}},
          color={95,95,95},
          thickness=0.5));
     end for;

        annotation (Line(points={{-37,50},{-26,50},{-26,49.8},{-0.2,49.8}}, color={0,0,127}));
    end EMBS_Body;

    model Nodes
      parameter Integer numNodes = 9;
      //parameter String SIDfileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
      EMBSlib.Types.Taylor origin[numNodes];
      EMBSlib.Types.Taylor phi[numNodes];
      EMBSlib.Types.Taylor psi[numNodes];
      EMBSlib.Types.Taylor AP[numNodes];

      parameter EMBSlib.SID_Data sid;

    algorithm



      for i in 1:numNodes loop
    /*
    //origin
    //===========================
  origin[i].order := EMBSlib.ExternalFunctions.getOriginOrder(sid, i);
  origin[i].nrow := EMBSlib.ExternalFunctions.getOriginNRows(sid, i);
  origin[i].ncol := EMBSlib.ExternalFunctions.getOriginNCols(sid, i);
  origin[i].nq := EMBSlib.ExternalFunctions.getOriginNq(sid, i);
  origin[i].nqn := EMBSlib.ExternalFunctions.getOriginNqn(sid, i);
  origin[i].structure := EMBSlib.ExternalFunctions.getOriginStructure(sid,
    i);
  origin[i].M0[1, 1] := EMBSlib.ExternalFunctions.getOriginM0(
    sid,
    i,
    1,
    1,
    origin[i].nrow,
    origin[i].ncol);
  origin[i].M0[2, 1] := EMBSlib.ExternalFunctions.getOriginM0(
    sid,
    i,
    2,
    1,
    origin[i].nrow,
    origin[i].ncol);
  origin[i].M0[3, 1] := EMBSlib.ExternalFunctions.getOriginM0(
    sid,
    i,
    3,
    1,
    origin[i].nrow,
    origin[i].ncol);
    for r0 in (1:origin[i].nrow) loop
      for q0 in (1:origin[i].nq) loop
        for c0 in (1:origin[i].ncol) loop
            origin[i].M1[r0,q0,c0] := EMBSlib.ExternalFunctions.getOriginM1(    sid,    i,    r0,q0,c0,    origin[i].nrow,    origin[i].nq,    origin[i].ncol);
    end for;
    end for;
    end for;

    //Phi
    //===========================
  phi[i].order := EMBSlib.ExternalFunctions.getPhiOrder(sid, i);
  phi[i].nrow := EMBSlib.ExternalFunctions.getPhiNRows(sid, i);
  phi[i].ncol := EMBSlib.ExternalFunctions.getPhiNCols(sid, i);
  phi[i].nq := EMBSlib.ExternalFunctions.getPhiNq(sid, i);
  phi[i].nqn := EMBSlib.ExternalFunctions.getPhiNqn(sid, i);
  phi[i].structure := EMBSlib.ExternalFunctions.getPhiStructure(sid,
    i);
  phi[i].M0[1, 1] := EMBSlib.ExternalFunctions.getPhiM0(
    sid,
    i,
    1,
    1,
    phi[i].nrow,
    phi[i].ncol);
  phi[i].M0[2, 1] := EMBSlib.ExternalFunctions.getPhiM0(
    sid,
    i,
    2,
    1,
    phi[i].nrow,
    phi[i].ncol);
  phi[i].M0[3, 1] := EMBSlib.ExternalFunctions.getPhiM0(
    sid,
    i,
    3,
    1,
    phi[i].nrow,
    phi[i].ncol);
    for r1 in (1:phi[i].nrow) loop
      for q1 in (1:phi[i].nq) loop
        for c1 in (1:phi[i].ncol) loop
            phi[i].M1[r1,q1,c1] := EMBSlib.ExternalFunctions.getPhiM1(    sid,    i,    r1,q1,c1,    phi[i].nrow,    phi[i].nq,    phi[i].ncol);
    end for;
    end for;
    end for;

    //Psi
    //===========================
  psi[i].order := EMBSlib.ExternalFunctions.getPsiOrder(sid, i);
  psi[i].nrow := EMBSlib.ExternalFunctions.getPsiNRows(sid, i);
  psi[i].ncol := EMBSlib.ExternalFunctions.getPsiNCols(sid, i);
  psi[i].nq := EMBSlib.ExternalFunctions.getPsiNq(sid, i);
  psi[i].nqn := EMBSlib.ExternalFunctions.getPsiNqn(sid, i);
  psi[i].structure := EMBSlib.ExternalFunctions.getPsiStructure(sid,
    i);
  psi[i].M0[1, 1] := EMBSlib.ExternalFunctions.getPsiM0(
    sid,
    i,
    1,
    1,
    psi[i].nrow,
    psi[i].ncol);
  psi[i].M0[2, 1] := EMBSlib.ExternalFunctions.getPsiM0(
    sid,
    i,
    2,
    1,
    psi[i].nrow,
    psi[i].ncol);
  psi[i].M0[3, 1] := EMBSlib.ExternalFunctions.getPsiM0(
    sid,
    i,
    3,
    1,
    psi[i].nrow,
    psi[i].ncol);
      for r2 in (1:psi[i].nrow) loop
      for q2 in (1:psi[i].nq) loop
        for c2 in (1:psi[i].ncol) loop
            psi[i].M1[r2,q2,c2] := EMBSlib.ExternalFunctions.getPsiM1(    sid,    i,    r2,q2,c2,    psi[i].nrow,    psi[i].nq,    psi[i].ncol);
    end for;
    end for;
    end for;

    //AP
    //===========================
  AP[i].order := EMBSlib.ExternalFunctions.getAPOrder(sid, i);
  AP[i].nrow := EMBSlib.ExternalFunctions.getAPNRows(sid, i);
  AP[i].ncol := EMBSlib.ExternalFunctions.getAPNCols(sid, i);
  AP[i].nq := EMBSlib.ExternalFunctions.getAPNq(sid, i);
  AP[i].nqn := EMBSlib.ExternalFunctions.getAPNqn(sid, i);
  AP[i].structure := EMBSlib.ExternalFunctions.getAPStructure(sid,
    i);
  AP[i].M0[1, 1] := EMBSlib.ExternalFunctions.getAPM0(
    sid,
    i,
    1,
    1,
    AP[i].nrow,
    AP[i].ncol);
  AP[i].M0[2, 1] := EMBSlib.ExternalFunctions.getAPM0(
    sid,
    i,
    2,
    1,
    AP[i].nrow,
    AP[i].ncol);
  AP[i].M0[3, 1] := EMBSlib.ExternalFunctions.getAPM0(
    sid,
    i,
    3,
    1,
    AP[i].nrow,
    AP[i].ncol);
    for r3 in (1:AP[i].nrow) loop
      for q3 in (1:AP[i].nq) loop
        for c3 in (1:AP[i].ncol) loop
            AP[i].M1[r3,q3,c3] := EMBSlib.ExternalFunctions.getAPM1(    sid,    i,    r3,q3,c3,    AP[i].nrow,    AP[i].nq,    AP[i].ncol);
    end for;
    end for;
    end for;
*/
      end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));
    end Nodes;

    model Node
      parameter EMBSlib.SID_Data sid;
      parameter Integer nq = 3;
      parameter Integer nodeArrayIdx = 1;
      parameter Integer nr0 = 3;


      //origin
      parameter Real origin_M0[nr0,1] = EMBSlib.ExternalFunctions.getM0Node( sid, "origin",nodeArrayIdx, nr0, 1) annotation (Evaluate=true);
      parameter Real origin_M1[nr0,nq,1] = EMBSlib.ExternalFunctions.getM1Node(sid, "origin",nodeArrayIdx, nr0, nq, 1) annotation (Evaluate=true);
      Real origin [nr0,1] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,1,origin_M0,origin_M1,q);

      //psi
      parameter Real psi_M0[nr0,nq] = EMBSlib.ExternalFunctions.getM0Node( sid, "psi",nodeArrayIdx, nr0, nq) annotation (Evaluate=true);
      parameter Real psi_M1[nr0,nq,nq] = EMBSlib.ExternalFunctions.getM1Node(sid, "psi",nodeArrayIdx, nr0, nq, nq) annotation (Evaluate=true);
      Real psi [nr0,nq] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,nq,psi_M0,psi_M1,q);

      //phi
      parameter Real phi_M0[nr0,nq] = EMBSlib.ExternalFunctions.getM0Node( sid, "phi",nodeArrayIdx, nr0, nq) annotation (Evaluate=true);
      parameter Real phi_M1[nr0,nq,nq] = EMBSlib.ExternalFunctions.getM1Node(sid, "phi",nodeArrayIdx, nr0, nq, nq) annotation (Evaluate=true);
      Real phi [nr0,nq] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,nq,phi_M0,phi_M1,q);

      Modelica.Blocks.Interfaces.RealInput q[nq] "modal coordinates"
        annotation (Placement(transformation(extent={{-122,-22},{-82,18}})));
      Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a
        annotation (Placement(transformation(extent={{-116,-72},{-84,-40}})));
      Modelica.SIunits.Force f[nr0] = frame_b.f "external force applied";

        Modelica.SIunits.Torque t[nr0] = frame_b.t "external torque applied";


      Modelica.SIunits.Force hde_i[nq] = transpose(phi)*f;
      Modelica.SIunits.Force hdt_i[nq] = f;


      parameter Modelica.SIunits.Diameter sphereDiameter=world.defaultBodyDiameter
        "Diameter of sphere" annotation (Dialog(
          tab="Animation",
          group="if animation = true",
          enable=animation));
                               Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape(r = frame_b.r_0- {1,0,0}*sphereDiameter/2,
     shapeType="sphere",
     length=0.1,
     width=0.1,
     height=0.1,
     lengthDirection={1,0,0})
     annotation (Placement(transformation(extent={{80,60},{100,80}})));



      Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b
        annotation (Placement(transformation(extent={{84,-74},{116,-42}})));
    protected
      outer Modelica.Mechanics.MultiBody.World world
        annotation (Placement(transformation(extent={{-98,60},{-78,80}})));
    equation
      Connections.branch(frame_a.R, frame_b.R);
      assert(cardinality(frame_a) > 0 or cardinality(frame_b) > 0,
        "Neither connector frame_a nor frame_b of FixedTranslation object is connected");

      frame_b.r_0 = frame_a.r_0 + origin[:,1];
      frame_b.R = frame_a.R;

      /* Force and torque balance */
      zeros(3) = frame_a.f + frame_b.f;
      zeros(3) = frame_a.t + frame_b.t;
      //connect(frame_a, frame_b) annotation (Line(      points={{-100,-56},{-2,-56},{-2,-58},{100,-58}},      color={95,95,95},      thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}})));
    end Node;
  end Components;

  class SID_Data
    extends ExternalObject;

    function constructor
      input String fileName;
      output SID_Data sid;
      external "C" sid = SID_Constructor(fileName)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});

    end constructor;

    function destructor
      input SID_Data sid;
      external "C" SID_Destructor(sid)
      annotation(Include="#include \"SID_Data.h\"",
                   Library={"readSIDlib"});
    end destructor;

  end SID_Data;

  package ExternalFunctions
    function getNumberOfNodes
      input EMBSlib.SID_Data sid;
      output Integer numNodes;
      external "C" numNodes = getNumberOfNodes(sid)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getNumberOfNodes;

    function getNumberOfModes
      input EMBSlib.SID_Data sid;
      output Integer numNodes;
      external "C" numNodes = getNumberOfModes(sid)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getNumberOfModes;

    //Origin Taylor============================================
    //=========================================================
    function getOriginOrder
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginOrder(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginOrder;

      function getOriginNRows
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numRows;
      external "C" numRows = getOriginNRows(sid, nodeIdx)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getOriginNRows;

    function getOriginNCols
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numCols;
    external"C" numCols = getOriginNCols(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginNCols;

    function getOriginNq
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginNq(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginNq;

      function getOriginNqn
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginNqn(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getOriginNqn;

     function getOriginStructure
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginStructure(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
     end getOriginStructure;

      function getM0ForNode
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer numRows;
      input Integer numCols;
      output Real[numRows, numCols] M0;
      external "C" getM0ArrforNode(sid, nodeIdx, M0)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getM0ForNode;

    function getOriginM0
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer c;
      input Integer dimR;
      input Integer dimC;
      output Real value;
      external "C" value = getOriginM0(sid,nodeIdx,r,c, dimR, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginM0;

      function getOriginM1
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer q;
      input Integer c;
      input Integer dimR;
      input Integer dimQ;
      input Integer dimC;
      output Real value;
      external "C" value = getOriginM1(sid,nodeIdx,r,q,c, dimR,dimQ, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getOriginM1;

    //Phi Taylor============================================
    //=========================================================
    function getPhiOrder
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPhiOrder(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPhiOrder;

      function getPhiNRows
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numRows;
      external "C" numRows = getPhiNRows(sid, nodeIdx)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getPhiNRows;

    function getPhiNCols
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numCols;
    external"C" numCols = getPhiNCols(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPhiNCols;

    function getPhiNq
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPhiNq(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPhiNq;

      function getPhiNqn
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPhiNqn(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getPhiNqn;

     function getPhiStructure
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPhiStructure(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
     end getPhiStructure;

    function getPhiM0
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer c;
      input Integer dimR;
      input Integer dimC;
      output Real value;
      external "C" value = getPhiM0(sid,nodeIdx,r,c, dimR, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPhiM0;

      function getPhiM1
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer q;
      input Integer c;
      input Integer dimR;
      input Integer dimQ;
      input Integer dimC;
      output Real value;
      external "C" value = getPhiM1(sid,nodeIdx,r,q,c, dimR,dimQ, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getPhiM1;

        //Psi Taylor============================================
    //=========================================================
    function getPsiOrder
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPsiOrder(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPsiOrder;

      function getPsiNRows
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numRows;
      external "C" numRows = getPsiNRows(sid, nodeIdx)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getPsiNRows;

    function getPsiNCols
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numCols;
    external"C" numCols = getPsiNCols(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPsiNCols;

    function getPsiNq
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPsiNq(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPsiNq;

      function getPsiNqn
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPsiNqn(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getPsiNqn;

     function getPsiStructure
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getPsiStructure(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
     end getPsiStructure;

    function getPsiM0
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer c;
      input Integer dimR;
      input Integer dimC;
      output Real value;
      external "C" value = getPsiM0(sid,nodeIdx,r,c, dimR, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getPsiM0;

      function getPsiM1
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer q;
      input Integer c;
      input Integer dimR;
      input Integer dimQ;
      input Integer dimC;
      output Real value;
      external "C" value = getPsiM1(sid,nodeIdx,r,q,c, dimR,dimQ, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getPsiM1;

    //AP Taylor================================================
    //=========================================================

    function getAPOrder
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getAPOrder(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getAPOrder;

      function getAPNRows
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numRows;
      external "C" numRows = getAPNRows(sid, nodeIdx)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getAPNRows;

    function getAPNCols
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numCols;
    external"C" numCols = getAPNCols(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getAPNCols;

    function getAPNq
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getAPNq(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getAPNq;

      function getAPNqn
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getAPNqn(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getAPNqn;

     function getAPStructure
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getAPStructure(sid, nodeIdx)
        annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
     end getAPStructure;

      function getAPM0
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer c;
      input Integer dimR;
      input Integer dimC;
      output Real value;
      external "C" value = getAPM0(sid,nodeIdx,r,c, dimR, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getAPM0;

      function getAPM1
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer q;
      input Integer c;
      input Integer dimR;
      input Integer dimQ;
      input Integer dimC;
      output Real value;
      external "C" value = getAPM1(sid,nodeIdx,r,q,c, dimR,dimQ, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getAPM1;

    function getMass
      input EMBSlib.SID_Data sid;
      output Real mass;
      external "C" mass = getMass(sid)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getMass;

    function getM0
      input EMBSlib.SID_Data sid;
      input String taylorName;
      input Integer nr;
      input Integer nc;
      output Real[nr,nc] m0;
      external "C" getM0(sid,taylorName,m0,nr,nc)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getM0;

    function getM1
      input EMBSlib.SID_Data sid;
      input String taylorName;
      input Integer nr;
      input Integer nq;
      input Integer nc;
      output Real[nr,nq,nc] m1;
      external "C" getM1(sid,taylorName,m1,nr,nq,nc)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getM1;

    function getM0Node
      input EMBSlib.SID_Data sid;
      input String taylorName;
      input Integer nodeIdx;
      input Integer nr;
      input Integer nc;
      output Real[nr,nc] m0;

      external "C" getM0Node(sid,taylorName,nodeIdx,m0,nr,nc)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getM0Node;

    function getM1Node
      input EMBSlib.SID_Data sid;
      input String taylorName;
      input Integer nodeIdx;
      input Integer nr;
      input Integer nq;
      input Integer nc;
      output Real[nr,nq,nc] m1;
      external "C" getM1Node(sid,taylorName,nodeIdx,m1,nr,nq,nc)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getM1Node;
  end ExternalFunctions;

  package Types
    record Taylor
      Integer order;
      Integer nrow;
      Integer ncol;
      Integer nq;
      Integer nqn;
      Integer structure;
      Real[3,3] M0; //nrows,ncol
      Real[3,3,3] M1; //rows,q,cols
    end Taylor;

    record Matrix31
        Real M[3];
    end Matrix31;
  end Types;

  package Functions
    function getOrigin
        input EMBSlib.SID_Data sid;
        input Integer nodeIdx;
        output EMBSlib.Types.Taylor nodeOrig;
    algorithm
      nodeOrig.order := EMBSlib.ExternalFunctions.getOriginOrder(sid, nodeIdx);
      nodeOrig.nrow := EMBSlib.ExternalFunctions.getOriginNRows(sid, nodeIdx);
      nodeOrig.ncol := EMBSlib.ExternalFunctions.getOriginNCols(sid, nodeIdx);
      nodeOrig.nq := EMBSlib.ExternalFunctions.getOriginNq(sid, nodeIdx);
      nodeOrig.nqn := EMBSlib.ExternalFunctions.getOriginNqn(sid, nodeIdx);
      nodeOrig.structure := EMBSlib.ExternalFunctions.getOriginStructure(sid,
        nodeIdx);
      nodeOrig.M0[1, 1] := EMBSlib.ExternalFunctions.getOriginM0(
            sid,
            nodeIdx,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.ncol);
      nodeOrig.M0[2, 1] := EMBSlib.ExternalFunctions.getOriginM0(
            sid,
            nodeIdx,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.ncol);
      nodeOrig.M0[3, 1] := EMBSlib.ExternalFunctions.getOriginM0(
            sid,
            nodeIdx,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.ncol);
      nodeOrig.M1[1, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            1,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[2, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            2,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[3, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            3,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[1, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            1,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[2, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            2,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[3, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            3,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[1, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            1,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[2, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            2,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[3, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            3,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
    end getOrigin;
  end Functions;

  model SimplePlate2 "Plate with force excitation"
    import FlexibleBodies;
    extends Modelica.Icons.Example;
    inner Modelica.Mechanics.MultiBody.World world(
      enableAnimation=true,
      n={0,0,-1},
      animateGravity=false,
      g=0,
      animateWorld=true)
      annotation (Placement(transformation(extent={{-112,10},{-92,30}},
            rotation=0)));

    Modelica.Mechanics.MultiBody.Forces.WorldForce worldForce(N_to_m=10)
      annotation (Placement(transformation(extent={{-2,-58},{18,-38}},
            rotation=0)));
    Modelica.Blocks.Sources.Sine sine(
      freqHz=4,
      phase=1.7,
      amplitude=0,
      offset=10)    annotation (Placement(transformation(extent={{-60,-50},{-40,
              -30}},     rotation=0)));
    Modelica.Blocks.Sources.Constant const(k=0)
      annotation (Placement(transformation(extent={{-60,-80},{-40,-60}},
            rotation=0)));
    FlexibleBodies.ModalBody ModalBody(
      wireFrameElasticScale=1,
      WavefrontFile=FlexibleBodies.DataDirectory + "cartopPragV3.obj",
      solidColor={155,155,155},
      wireFrameColor={155,0,0},
      nodeAxesAnimation=true,
      SID_fileName=FlexibleBodies.DataDirectory + "cartopPragV32.SID_FEM",
      solidElasticScale=1,
      Nodes={1,2,12,26,78,84,107,143,149},
      solidAnimation=false,
      wireFrameAnimation=true,
      fontSize=8,
      showNodeNumbers=true,
      get_c=true,
      get_u=true,
      get_phi=true,
      native=true) annotation (Placement(transformation(extent={{-20,-10},{40,50}},
            rotation=0)));


      //reverse engineering
      parameter Integer nq = 3;
      Real q[nq] = ModalBody.q;
      Real q_d[nq] = der(q);
      Real q_dd[nq] = der(q_d);
      Real m = ModalBody.nBody.modal.mass;
      Real mI[3,3] = identity(3)*m;

      //mdCM
      Real mdCM_M0[3,1] = ModalBody.nBody.modal.mdCM.M0;
      Real mdCM_M1[3,nq,1] = ModalBody.nBody.modal.mdCM.M1;
      Real mdCM_ [3,1] = EMBSlib.MatrixFunctions.getTaylorFunction(3,nq,1,mdCM_M0,mdCM_M1,q);
      //Real mdCM [3,3] = {{mdCM_[1,1],0,0},{0,mdCM_[2,1],0},{0,0,mdCM_[3,1]}};
      Real cm[3,1] = mdCM_/m;
      Real mdCM_tilde[3,3] = {{0, -mdCM_[3,1],mdCM_[2,1]},  {mdCM_[3,1],0,-mdCM_[1,1]}, {-mdCM_[2,1],mdCM_[1,1],0}};

      //J
      Real J_M0[6,1] = ModalBody.nBody.modal.J.M0;
      Real J_M1[6,nq,1] = ModalBody.nBody.modal.J.M1;
      Real J_ [6,1] = EMBSlib.MatrixFunctions.getTaylorFunction(6,nq,1,J_M0,J_M1,q);
      Real J [3,3] = {{J_[1,1],J_[4,1],J_[5,1]},{J_[4,1],J_[2,1],J_[6,1]},{J_[5,1],J_[6,1],J_[3,1]}};
      //Real J[6,1] =  J_M1*q;

      //Ct
      Real Ct_M0[3,3] = ModalBody.nBody.modal.Ct.M0;
      Real Ct_M1[3,nq,3] = ModalBody.nBody.modal.Ct.M1;
      Real Ct [3,3] = EMBSlib.MatrixFunctions.getTaylorFunction(3,nq,3,Ct_M0,Ct_M1,q);

      //Cr
      Real Cr_M0[3,3] = ModalBody.nBody.modal.Cr.M0;
      Real Cr_M1[3,nq,3] = ModalBody.nBody.modal.Cr.M1;
      Real Cr [3,3] = EMBSlib.MatrixFunctions.getTaylorFunction(3,nq,3,Cr_M0,Cr_M1,q);

      //Gr
      Real Gr_ [3,9] = ModalBody.nBody.modal.Gr;
      Real Gr [3,3] = EMBSlib.MatrixFunctions.getGrMatrix(3,nq,9,Gr_,q_d);//! q_d

      //Ge
      Real Ge_ [3,9] = ModalBody.nBody.modal.Ge;
      Real Ge [3,3] = EMBSlib.MatrixFunctions.getGrMatrix(3,nq,9,Ge_,q_d);//! q_d

      //Me
      Real Me[3,3] = identity(3);

      //Oe
      Real Oe_M0[3,6] = ModalBody.nBody.modal.Oe.M0;
      Real Oe_M1[3,nq,6] = ModalBody.nBody.modal.Oe.M1;
      Real Oe [nq,6] = EMBSlib.MatrixFunctions.getTaylorFunction(3,nq,6,Oe_M0,Oe_M1,q);

      //OMega
      Real OMega[6] = {phi_d[1]^2, phi_d[2]^2, phi_d[3]^2, phi_d[1]*phi_d[2], phi_d[2]*phi_d[3], phi_d[1]*phi_d[3]};

      //ksigma
      Real ksigma[3,1] = zeros(3,1);

      //Ke
      Real Ke[nq,nq] = ModalBody.nBody.modal.Ke;

      //De
      Real De[nq,nq] =  ModalBody.nBody.modal.De;

      // references
      Real f[3] = ModalBody.frame_ref.f;
      Real t[3] = ModalBody.frame_ref.t;
      Real t_diff[3] = -ModalBody.frame_ref.t + h_r + t_comp; // this should be zero

      //external loads on all nodes (without frame of reference)
      Real h_t[3] = sum(ModalBody.nodes[i].f for i in 1:ModalBody.n_Nodes);
      Real h_r[3] = sum(cross(ModalBody.nodes[i].f,ModalBody.nodes[i].r_0) for i in 1:ModalBody.n_Nodes);

      //kinematic equations
      Real r_0[3] = ModalBody.frame_ref.r_0;
      Real r_0_d[3] = der(r_0);
      Real r_0_dd[3] = der(r_0_d);
      Real phi_d[3] = Modelica.Mechanics.MultiBody.Frames.angularVelocity2(ModalBody.frame_ref.R);
      Real phi_dd[3] = der(phi_d);
      Real omega_tilde[3,3] = {{0, -phi_d[3],phi_d[2]},  {phi_d[3],0,-phi_d[1]}, {-phi_d[2],phi_d[1],0}};

      // kinematic equations --> translation
      Real M_t[3] = mI*r_0_dd + transpose(mdCM_tilde)*phi_dd + transpose(Ct)*q_dd;
      Real k_omega_t[3] = 2*omega_tilde*transpose(Ct)*q_d + res2_1[:,1];
      Real res2_1[3,1] = omega_tilde*omega_tilde*cm; // need the intermediate value to fix the dimensions
      Real f_comp[3] = M_t+k_omega_t;
      Real hd_t[3] = sum(identity(3)*ModalBody.nodes[i].f for i in 1:ModalBody.n_Nodes);
      Real residual_t[3] = f_comp - hd_t - ModalBody.frame_ref.f;  // this should be zero

      // kinematic equations --> rotation
      Real M_r[3] = mdCM_tilde*r_0_dd + J*phi_dd + transpose(Cr) *q_dd;
      Real k_omega_r[3] = Gr*phi_dd + omega_tilde*J*phi_d;
      Real t_comp[3] = M_r+k_omega_r;
      Real hd_r[3] = sum(identity(3)*ModalBody.nodes[i].f for i in 1:ModalBody.n_Nodes);

      //kinematic equations --> modal
      Real M_q[3] = Ct*r_0_dd + Cr*phi_dd + Me*q_dd;
      Real k_omega_q[3] =  Ge*phi_d + Oe*OMega;
      Real k_q[nq] = ksigma[:,1] + Ke*q +De*q_d;
      Real q_comp[nq] = M_q + k_omega_q + k_q;
      Real hde[nq] = phi*ModalBody.nodes[forceNode].f;

      //node position

      //Verschiebung
      Real u[3]= ModalBody.nBody.u_out[forceNode, :];//forceNode
      //phi of force node
      Real phi_M0[3,3] = ModalBody.nBody.modal.nodes[forceNode].phi.M0;
      Real phi_M1[3,nq,3] = ModalBody.nBody.modal.nodes[forceNode].phi.M1;
      Real phi [3,3] = EMBSlib.MatrixFunctions.getTaylorFunction(3,nq,3,phi_M0,phi_M1,q);
      //psi of force node
      Real psi_M0[3,3] = ModalBody.nBody.modal.nodes[forceNode].psi.M0;
      Real psi_M1[3,nq,3] = ModalBody.nBody.modal.nodes[forceNode].psi.M1;
      Real psi [3,3] = EMBSlib.MatrixFunctions.getTaylorFunction(3,nq,3,psi_M0,psi_M1,q);


      Real uForceNode[3] = phi*q;

      parameter Integer forceNode = 9;//5 works?!


  algorithm

  equation


    connect(sine.y, worldForce.force[3]) annotation (Line(points={{-39,-40},{
            -20,-40},{-20,-46.6667},{-4,-46.6667}},  color={0,0,127}));
    connect(const.y, worldForce.force[1]) annotation (Line(points={{-39,-70},{
            -20,-70},{-20,-49.3333},{-4,-49.3333}},   color={0,0,127}));
    connect(const.y, worldForce.force[2]) annotation (Line(points={{-39,-70},{-20,
            -70},{-20,-48},{-4,-48}},       color={0,0,127}));
    connect(worldForce.frame_b, ModalBody.nodes[forceNode]) annotation (Line(
        points={{18,-48},{60,-48},{60,20},{40,20}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    connect(world.frame_b, ModalBody.frame_ref) annotation (Line(
        points={{-92,20},{-20,20}},
        color={95,95,95},
        thickness=0.5));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}})),
                         Documentation(info="<html>
<p>
The flexible structure of this academic example is a rectangular plate. Its displacement field is approximated with three eigenmodes in the frequency range below 3 Hz. Three corners of the plate are simply supported, free boundary conditions have been applied to the fourth corner. These properties are the result of the eigenvalue analysis, that have been performed with the finite element plate model in Ansys with the corresponding displacement constraints imposed. The origin of the reference frame coincides with one corner of the plate.
</p>
 
<p>
A force is applied at node 149 normal to the plate and exites the linear structure harmonically. The picture below shows the deformation field twice: the grey structure depicts the in scale displacments, the red wireframe exaggerates the elastic displacements by a factor of 30. Each clamped node that is specified is visualised  by a small, green coordinate system.
</p>
<br clear=\"all\">
<div align=\"center\">
<img
 align=\"bottom\" border=\"0\"
 src=\"modelica://FlexibleBodies/Resources/Images/Examples/PlateAnimation.png\" width=\"533\"
 ><br>
Animation model SimplePlate
</div>
<p>
The plot below shows the transient modal amplitudes of the three modes due to the force excitation.
</p>

<br clear=\"all\">
<div align=\"center\">
<img
 align=\"bottom\" border=\"0\"
 src=\"modelica://FlexibleBodies/Resources/Images/Examples/PlatePlot.png\"
 ><br>
Results of the model SimplePlate
</div>
</html>",
  revisions="<html>
<table border=0 cellspacing=0 cellpadding=2>
  <tr>
    <td valign=\"center\"><img src=\"modelica://FlexibleBodies/Resources/Images/dlr_logo.png\"  width=60 ></td>
    <td valign=\"center\">&nbsp;&nbsp;&nbsp;</td>
    <td valign=\"center\"><b>Copyright</b>
      <br><b>&copy; 1999-2012, DLR Institute of Robotics and Mechatronics</b>
      <br><b>&copy; 2012-2019, DLR Institute of System Dynamics and Control</b></td>
  </tr>
</table>
</html>"),   experiment(StopTime=1.01, Tolerance=1e-006),
      experimentSetupOutput,
      Commands(file="modelica://FlexibleBodies/Resources/Scripts/Settings.mos"
          "Settings",
        file="modelica://FlexibleBodies/Resources/Scripts/SimulateandPlotSimplePlate.mos"
          "Simulate and Plot Results",
        file="modelica://FlexibleBodies/Resources/Scripts/PlotSimplePlate.mos"
          "Plot Results"));
  end SimplePlate2;

  package MatrixFunctions
    function getTaylorFunction
          input Integer nr;
          input Integer nq;
          input Integer nc;
          input Real[nr,nc] M0;
          input Real[nr,nq,nc] M1;
          input Real[nq] q;
          output Real[nr,nc] M;
    protected
          Real[nr,nc] aux;
    algorithm
          M := zeros(nr,nc);
          for i in 1:nq loop
            aux := M1[:,i,:]*q[i];
            M := M+aux;
          end for;
          M := M0 +M;
    end getTaylorFunction;

        function getGrMatrix
          input Integer nr;//(for Gr=3 for Ge=nq)
          input Integer nq;//nq
          input Integer nc;//3*nq
          input Real[nr,nc] Min;
          input Real[nq] q;
          output Real[nq, 3] Mout;
    protected
          Real[3,3] aux;
          Integer c=1;
        algorithm
          Mout := zeros(3,3);
          for i in 1:nq loop
            c :=(i - 1)*3;//column index
            aux[:,1] := Min[:,c+1];
            aux[:,2] := Min[:,c+2];
            aux[:,3] := Min[:,c+3];
            Mout := Mout+aux*q[i];
          end for;
        end getGrMatrix;
  end MatrixFunctions;
  annotation (uses(Modelica(version="3.2.3")));
end EMBSlib;
