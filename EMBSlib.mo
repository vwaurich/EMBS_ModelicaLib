within ;
package EMBSlib

  model Balken
    extends Modelica.Icons.Example;

    Components.EMBS_Body eMBS_Body(
      numNodes=6,
      numModes=11,
      SIDfileName=Modelica.Utilities.Files.loadResource(
          "modelica://EMBSlib/Resources/Data/Balken.SID_FEM"))
        annotation (Placement(transformation(extent={{-20,0},{0,20}})));
    inner Modelica.Mechanics.MultiBody.World world(g=1,
      n(displayUnit="1") = {0,0,-1})
        annotation (Placement(transformation(extent={{-100,0},{-80,20}})));

    Modelica.Mechanics.MultiBody.Forces.WorldForce force
      annotation (Placement(transformation(extent={{6,-48},{26,-28}})));
    Modelica.Blocks.Sources.RealExpression fIn[3](y={0,1e5*time,0})
      annotation (Placement(transformation(extent={{-86,-58},{-66,-38}})));
  equation
    connect(eMBS_Body.frame_ref, world.frame_b) annotation (Line(
        points={{-20,10},{-80,10}},
        color={95,95,95},
        thickness=0.5));

    connect(fIn.y, force.force) annotation (Line(points={{-65,-48},{-32,-48},{-32,
            -38},{4,-38}}, color={0,0,127}));
    connect(force.frame_b, eMBS_Body.frame_node[2]) annotation (Line(
        points={{26,-38},{36,-38},{36,-40},{52,-40},{52,9.2},{0,9.2}},
        color={95,95,95},
        thickness=0.5));
    annotation (experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
  end Balken;

  model EMBS_bodyExample
    extends Modelica.Icons.Example;

    Components.EMBS_Body eMBS_Body(SIDfileName=
          Modelica.Utilities.Files.loadResource(
          "modelica://EMBSlib/Resources/Data/cartopPragV32.SID_FEM"))
        annotation (Placement(transformation(extent={{-20,0},{0,20}})));

    inner Modelica.Mechanics.MultiBody.World world(g=9.81,
      n(displayUnit="1") = {0,0,-1})
        annotation (Placement(transformation(extent={{-100,0},{-80,20}})));

    Modelica.Mechanics.MultiBody.Joints.Revolute revFix(n(displayUnit="1") = {1,
        0,0}, w(fixed=true, start=0))
      annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
  equation
    connect(world.frame_b, revFix.frame_a) annotation (Line(
        points={{-80,10},{-58,10}},
        color={95,95,95},
        thickness=0.5));
    connect(eMBS_Body.frame_ref, revFix.frame_b) annotation (Line(
        points={{-20,10},{-38,10}},
        color={95,95,95},
        thickness=0.5));
    annotation (experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
  end EMBS_bodyExample;

  model TestSIDRead
        parameter Integer numNodes = 9;
        parameter Integer numModes = 3 "the number of modes given in the SID file";
        parameter String SIDfileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
        constant Integer nr0 = 3;
        parameter Integer nq = numModes;

        parameter EMBSlib.SID_File sid = EMBSlib.SID_File(SIDfileName) annotation (Evaluate=true);
        //mass
        parameter Modelica.SIunits.Mass mass = EMBSlib.ExternalFunctions_C.getMass(sid) annotation(Evaluate=true);
        parameter Modelica.SIunits.Mass mI[nr0,nr0] = identity(nr0)*mass;

        //mdCM
        parameter Real mdCM_M0[nr0,1] = EMBSlib.ExternalFunctions_C.getM0( sid, "mdCM", nr0, 1) annotation (Evaluate=true);
        parameter Real mdCM_M1[nr0,nq,1] = EMBSlib.ExternalFunctions_C.getM1(sid,   "mdCM", nr0, nq, 1) annotation (Evaluate=true);
        Real mdCM [nr0,1] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,1,mdCM_M0,mdCM_M1,q);

        Real q[nq] = {0.1,0.1,0.1};



    Components.Node node(sid=sid)
      annotation (Placement(transformation(extent={{32,0},{52,22}})));
    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-48,-24},{-28,-4}})));
    Modelica.Blocks.Sources.RealExpression qIn[3](y=q)
      annotation (Placement(transformation(extent={{-66,42},{-46,62}})));
  equation
    connect(world.frame_b, node.frame_a) annotation (Line(
        points={{-28,-14},{32,-14},{32,4.84}},
        color={95,95,95},
        thickness=0.5));
    connect(qIn.y, node.q) annotation (Line(points={{-45,52},{-6,52},{-6,10.78},{31.8,
            10.78}}, color={0,0,127}));
          annotation (Line(points={{-37,50},{-26,50},{-26,49.8},{-0.2,49.8}}, color={0,0,127}));
  end TestSIDRead;

  package Components
    model EMBS_Body
          parameter Integer numNodes = 9;
          parameter Integer numModes = 3 "the number of modes given in the SID file";
          parameter String SIDfileName =  Modelica.Utilities.Files.loadResource("modelica://EMBSlib/Resources/Data/cartopPragV32.SID_FEM");

          parameter Boolean animation = true annotation (Dialog(    tab="Animation"));
          parameter Real deformationScalingFactor= 1 annotation(Dialog(tab="Animation"));
          parameter Real sphereDiameter = 0.1 annotation(Dialog(tab="Animation"));
          parameter Real coordinateSystemScalingFactor = 0.2 annotation(Dialog(tab="Animation"));

          constant Integer nr0 = 3 "dimension in space";

          Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_ref
         annotation (Placement(transformation(extent={{-116,-16},{-84,16}})));
          Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_node[numNodes] annotation (Placement(transformation(extent={{84,-16},
                    {116,16}})));

          Node nodes[numNodes](each sid=sid, each nq=nq, nodeArrayIdx = 1:numNodes, each deformationScalingFactor=deformationScalingFactor, each sphereDiameter=sphereDiameter) annotation (Placement(transformation(extent={{-10,-4},{10,16}})));

          Real q[nq]( start=zeros(nq))  "modal coordinates";
          Real qd[nq]( start=zeros(nq));
          Real qdd[nq] = der(qd);

          Modelica.SIunits.Position r_0[3](start={0,0,0}) = frame_ref.r_0  "position of frame of reference";
          Modelica.SIunits.Velocity v[3] = der(r_0) "velocity of frame of reference";
          Modelica.SIunits.Acceleration a[3] = der(v) "acceleration of frame of reference";
          Modelica.SIunits.Acceleration g_0[3] =  Modelica.Mechanics.MultiBody.Frames.resolve2(frame_ref.R,world.gravityAcceleration(frame_ref.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_ref.R,
        cm[:,1]))) "Gravity acceleration resolved in world frame" annotation(Evaluate=true);

          Modelica.SIunits.AngularVelocity omega[3] = Modelica.Mechanics.MultiBody.Frames.angularVelocity2(frame_ref.R);
          Modelica.SIunits.AngularAcceleration omega_d[3] = der(omega);
          Modelica.SIunits.AngularVelocity omega_tilde[3,3] = {{0, -omega[3],omega[2]},  {omega[3],0,-omega[1]}, {-omega[2],omega[1],0}};

          // kinematic equations --> translation
          Real M_t[3] = mI*(a-g_0) + transpose(mdCM_tilde)* omega_d + transpose(Ct)* qdd;
          Real k_omega_t[3] = 2*omega_tilde*transpose(Ct)*qd + res2_1[:,1];
          Real res2_1[3,1] = omega_tilde*omega_tilde*mdCM; // need the intermediate value to fix the dimensions
          Real hd_t[3] = sum(identity(3)*nodes[i].f for i in 1:numNodes);

          // kinematic equations --> rotation
          Real M_r[3] = mdCM_tilde*(a-g_0) + J* omega_d + transpose(Cr) * qdd;
          Real k_omega_r[3] = Gr* omega + omega_tilde*J*omega;
          Real k_omega_r_1[3] = Gr* omega;
          Real k_omega_r_2[3] = omega_tilde*J*omega;

          //Real t_comp[3] = M_r+k_omega_r;
          Real hd_r[3] = sum(identity(3)*nodes[i].t for i in 1:numNodes);

          //kinematic equations --> modal
          Real M_q[nq] = Ct*(a-g_0) + Cr* omega_d + Me* qdd;
          Real k_omega_q[nq] =  Ge*omega + Oe_*OMega;
          Real k_omega_q_1[nq] =  Ge*omega;
          Real k_omega_q_2[nq] =  Oe_*OMega;

          Real k_q[nq] = ksigma[:,1] + Ke*q +De*qd;
          Real hd_e[nq] = sum(nodes[i].hde_i for i in 1:numNodes);

          Modelica.Blocks.Sources.RealExpression qExp[nq](y=q)
            annotation (Placement(transformation(extent={{-58,40},{-38,60}})));

          Modelica.SIunits.Torque[3] t_rest;
          Modelica.SIunits.Force[3] f_rest;
          Modelica.Mechanics.MultiBody.Forces.WorldForce force(resolveInFrame=
            Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.frame_b, N_to_m=0.1)
        annotation (Placement(transformation(extent={{-86,18},{-66,38}})));
          Modelica.Blocks.Sources.RealExpression f_elast1[nr0](y=f_rest)
        annotation (Placement(transformation(extent={{-124,16},{-104,36}})));
          Modelica.Mechanics.MultiBody.Forces.WorldTorque torque(resolveInFrame=
            Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.frame_b, Nm_to_m=0.1)
        annotation (Placement(transformation(extent={{-88,46},{-68,66}})));
          Modelica.Blocks.Sources.RealExpression t_elast2[nr0](y=t_rest)
        annotation (Placement(transformation(extent={{-124,46},{-104,66}})));
          Modelica.Mechanics.MultiBody.Visualizers.FixedFrame fixedFrame[numNodes](each length=
            coordinateSystemScalingFactor)
        annotation (Placement(transformation(extent={{80,38},{100,58}})));

          //SID File Data

    protected
          parameter EMBSlib.SID_File sid = EMBSlib.SID_File(SIDfileName) annotation (Evaluate=true);
          //mass
          parameter Modelica.SIunits.Mass mass = EMBSlib.ExternalFunctions_C.getMass(sid) annotation(Evaluate=true);
          parameter Modelica.SIunits.Mass mI[nr0,nr0] = identity(nr0)*mass;

          //mdCM
          parameter Real mdCM_M0[nr0,1] = EMBSlib.ExternalFunctions_C.getM0( sid, "mdCM", nr0, 1) annotation (Evaluate=true);
          parameter Real mdCM_M1[nr0,nq,1] = EMBSlib.ExternalFunctions_C.getM1(sid,   "mdCM", nr0, nq, 1) annotation (Evaluate=true);
          Real mdCM [nr0,1] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,1,mdCM_M0,mdCM_M1,q);
          Real cm[nr0,1] = mdCM/mass "centre of mass";
          Real mdCM_tilde[nr0,nr0] = {{0, -mdCM[3,1],mdCM[2,1]},  {mdCM[3,1],0,-mdCM[1,1]}, {-mdCM[2,1],mdCM[1,1],0}};

          //J, is always in a special structure in [6,1]
          constant Integer nJ = 6 "J is always a vektor of size 6";
          parameter Real J_M0[nJ,1]= EMBSlib.ExternalFunctions_C.getM0( sid, "J", nJ, 1) annotation (Evaluate=true);
          parameter Real J_M1[nJ,nq,1]= EMBSlib.ExternalFunctions_C.getM1( sid, "J", nJ, nq, 1) annotation (Evaluate=true);
          Real J_ [nJ,1] = EMBSlib.MatrixFunctions.getTaylorFunction(nJ,nq,1,J_M0,J_M1,q);
          Real J [nr0,nr0] = {{J_[1,1],J_[4,1],J_[5,1]},{J_[4,1],J_[2,1],J_[6,1]},{J_[5,1],J_[6,1],J_[3,1]}};

          //Ct
          parameter Real Ct_M0[nq,nr0] = EMBSlib.ExternalFunctions_C.getM0( sid, "Ct", nq, nr0) annotation (Evaluate=true);
          parameter Real Ct_M1[nq,nq,nr0] = EMBSlib.ExternalFunctions_C.getM1( sid, "Ct", nq, nq, nr0) annotation (Evaluate=true);
          Real Ct [nq,nr0] = EMBSlib.MatrixFunctions.getTaylorFunction(nq,nq,nr0,Ct_M0,Ct_M1,q);

          //Cr
          parameter Real Cr_M0[nq,nr0] = EMBSlib.ExternalFunctions_C.getM0( sid, "Cr", nq, nr0) annotation (Evaluate=true);
          parameter Real Cr_M1[nq,nq,nr0] = EMBSlib.ExternalFunctions_C.getM1( sid, "Cr", nq,  nq, nr0) annotation (Evaluate=true);
          Real Cr [nq,nr0] = EMBSlib.MatrixFunctions.getTaylorFunction(nq,nq,nr0,Cr_M0,Cr_M1,q);

          //Gr
          parameter Real Gr_ [nr0,nr0*nq] = EMBSlib.ExternalFunctions_C.getM0( sid, "Gr", nr0, nr0*nq) annotation (Evaluate=true);
          Real Gr [nr0,nr0] = EMBSlib.MatrixFunctions.getGrMatrix(nr0,nq,nr0*nq,Gr_,qd);//! this is the sum of the product of the GR matrices with the respective qd_i

          //Ge
          parameter Real Ge_ [nq,nr0*nq] = EMBSlib.ExternalFunctions_C.getM0( sid, "Ge", nq, nr0*nq) annotation (Evaluate=true);
          Real Ge [nq,nr0] = EMBSlib.MatrixFunctions.getGeMatrix(nq,nq,nr0*nq,Ge_,qd);//! this is the sum of the product of the Ge matrices with the respective qd_i

          //Me
          parameter Real Me[nq,nq] = identity(nq) annotation (Evaluate=true);

          //Oe
          constant Integer nOe = 6 "its always 6";
          parameter Real Oe_M0[nq,nOe] = EMBSlib.ExternalFunctions_C.getM0( sid, "Oe", nq, nOe) annotation (Evaluate=true);
          parameter Real Oe_M1[nq,nq,nOe] = EMBSlib.ExternalFunctions_C.getM1( sid, "Oe", nq,  nq, nOe) annotation (Evaluate=true);
          Real Oe_[nq,nOe] = EMBSlib.MatrixFunctions.getTaylorFunction(nq,nq,6,Oe_M0,Oe_M1,q);//wallrapp says on page 12 table 1 there is acompact form, not used here?!?!?!?

          //OMega
          Real OMega[nOe] = {omega[1]^2, omega[2]^2, omega[3]^2, omega[1]*omega[2], omega[2]*omega[3], omega[1]*omega[3]};

          //ksigma
          parameter Real ksigma[nq,1] = zeros(nq,1) annotation (Evaluate=true); //!! For simplicity....change this

          //Ke
          parameter Real Ke[nq,nq] = EMBSlib.ExternalFunctions_C.getM0( sid, "Ke", nq, nq) annotation (Evaluate=true);

           //De
          parameter Real De[nq,nq] =  EMBSlib.ExternalFunctions_C.getM0( sid, "De", nq, nq) annotation (Evaluate=true);

          parameter Integer nq = numModes;

          outer Modelica.Mechanics.MultiBody.World world
        annotation (Placement(transformation(extent={{-100,80},{-80,100}})));

    equation
         qd = der(q);

         //kinematic equations
         M_t +k_omega_t = -f_rest;
         M_r + k_omega_r = -t_rest;
         M_q + k_omega_q + k_q = hd_e;

         for i in 1:numNodes loop
          connect(qExp.y, nodes[i].q)    annotation (Line(points={{-37,50},{-18,50},{-18,5.8},{-10.2,5.8}},  color={0,0,127}));
          connect(nodes[i].frame_b, frame_node[i]) annotation (Line(
              points={{10,0.2},{25,0.2},{25,0},{100,0}},
              color={95,95,95},
              thickness=0.5));
          connect(nodes[i].frame_a, frame_ref) annotation (Line(
          points={{-10,0.4},{-55,0.4},{-55,0},{-100,0}},
          color={95,95,95},
          thickness=0.5));
         end for;

      connect(force.frame_b, frame_ref) annotation (Line(
          points={{-66,28},{-64,28},{-64,0},{-100,0}},
          color={95,95,95},
          thickness=0.5));
      connect(f_elast1.y, force.force) annotation (Line(points={{-103,26},{-96,
              26},{-96,28},{-88,28}}, color={0,0,127}));
      connect(torque.frame_b, force.frame_b) annotation (Line(
          points={{-68,56},{-68,28},{-66,28}},
          color={95,95,95},
          thickness=0.5));
      connect(t_elast2.y, torque.torque)
        annotation (Line(points={{-103,56},{-90,56}}, color={0,0,127}));

      connect(fixedFrame.frame_a, frame_node) annotation (Line(
          points={{80,48},{62,48},{62,0},{100,0}},
          color={95,95,95},
          thickness=0.5));
            annotation (Line(points={{-37,50},{-26,50},{-26,49.8},{-0.2,49.8}}, color={0,0,127}), Icon(
            graphics={
            Line(points={{-96,0},{-58,44},{-24,0},{-6,52},{-58,44},{-36,82},{-6,
                  52},{24,60},{-36,82},{52,78},{24,60},{56,20},{-6,52},{24,-18},
                  {-24,0},{-96,0},{-60,-42},{-24,0},{-20,-32},{-60,-42},{-32,
                  -66},{-20,-32},{8,-54},{-32,-66},{44,-76},{8,-54},{24,-18},{
                  44,-76},{68,-32},{24,-18},{56,20},{52,78},{86,28},{56,20},{68,
                  -32},{98,0},{86,28},{68,-32}}, color={28,108,200}),
            Line(points={{-20,-32},{24,-18}}, color={28,108,200}),
            Text(
              extent={{-64,118},{68,84}},
              lineColor={28,108,200},
              textString="%numNodes nodes"),
            Text(
              extent={{-70,-74},{78,-116}},
              lineColor={28,108,200},
              textString="%numModes modes")}));
    end EMBS_Body;

    model Node
          parameter EMBSlib.SID_File sid;
          parameter Integer nq = 3;
          parameter Integer nodeArrayIdx = 1;
          parameter Integer nr0 = 3;

          parameter Real deformationScalingFactor= 1;
          parameter Boolean animation = true annotation (Dialog(
              tab="Animation"));

          //origin (if bodies are attached to the nodes, origin has to be derived due to index reduction, keeping it as a parameter works)
          parameter Real origin_M0[nr0,1] = EMBSlib.ExternalFunctions_C.getM0Node( sid, "origin",nodeArrayIdx, nr0, 1) annotation (Evaluate=true);
          parameter Real origin_M1[nr0,nq,1] = EMBSlib.ExternalFunctions_C.getM1Node(sid, "origin",nodeArrayIdx, nr0, nq, 1) annotation (Evaluate=true);
          Real origin_[nr0,1] = origin_M0;//EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,1,origin_M0,origin_M1,q);
          Real origin[nr0] = origin_[:,1];

          //psi
          parameter Real psi_M0[nr0,nq] = EMBSlib.ExternalFunctions_C.getM0Node( sid, "psi",nodeArrayIdx, nr0, nq) annotation (Evaluate=true);
          parameter Real psi_M1[nr0,nq,nq] = EMBSlib.ExternalFunctions_C.getM1Node(sid, "psi",nodeArrayIdx, nr0, nq, nq) annotation (Evaluate=true);
          Real psi [nr0,nq] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,nq,psi_M0,psi_M1,q);
          Modelica.SIunits.Angle theta [nr0] = psi*q "elastic rotation";

          //phi(if bodies are attached to the nodes, phi has to be derived due to index reduction, keeping it as a parameter works)
          parameter Real phi_M0[nr0,nq] = EMBSlib.ExternalFunctions_C.getM0Node( sid, "phi",nodeArrayIdx, nr0, nq) annotation (Evaluate=true);
          parameter Real phi_M1[nr0,nq,nq] = EMBSlib.ExternalFunctions_C.getM1Node(sid, "phi",nodeArrayIdx, nr0, nq, nq) annotation (Evaluate=true);
          parameter Real phi [nr0,nq] = phi_M0;//EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,nq,phi_M0,phi_M1,q);
          Modelica.SIunits.Position u[nr0] = phi*q "elastic displacement" annotation(each stateSelect=StateSelect.never);


          //AP
          parameter Real AP_M0[nr0,nr0] = EMBSlib.ExternalFunctions_C.getM0Node( sid, "AP", nodeArrayIdx, nr0, nr0) annotation (Evaluate=true);
          parameter Real AP_M1[nr0,nq,nr0] = EMBSlib.ExternalFunctions_C.getM1Node(sid, "AP",nodeArrayIdx, nr0, nq, nr0) annotation (Evaluate=true);
          Real AP [nr0,nr0] = EMBSlib.MatrixFunctions.getTaylorFunction(nr0,nq,nr0,AP_M0,zeros(nr0,nq,nr0),q);

          Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a
            annotation (Placement(transformation(extent={{-116,-72},{-84,-40}})));

          Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b
            annotation (Placement(transformation(extent={{84,-74},{116,-42}})));

          Modelica.Blocks.Interfaces.RealInput q[nq] "modal coordinates"
            annotation (Placement(transformation(extent={{-122,-22},{-82,18}})));
          Real q_d[nq] = der(q);

          Modelica.SIunits.Force f[nr0] = frame_b.f "external force applied";

          Modelica.SIunits.Torque t[nr0] = frame_b.t "external torque applied";

          Modelica.SIunits.Force hde_i[nq] = transpose(phi)*f+transpose(psi)*t;

          parameter Modelica.SIunits.Diameter sphereDiameter=world.defaultBodyDiameter
            "Diameter of sphere" annotation (Dialog(
              tab="Animation",
              group="if animation = true",
              enable=animation));

          Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape(
        r=frame_b.r_0,
         shapeType="sphere",
         length=sphereDiameter,
         width=sphereDiameter,
         height=sphereDiameter,
         lengthDirection={1,0,0})
         annotation (Placement(transformation(extent={{80,60},{100,80}})));




    protected
          outer Modelica.Mechanics.MultiBody.World world
            annotation (Placement(transformation(extent={{-98,60},{-78,80}})));
          Modelica.Mechanics.MultiBody.Frames.Orientation R_theta = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3},theta,zeros(3));
    equation
          Connections.branch(frame_a.R, frame_b.R);
          assert(cardinality(frame_a) > 0 or cardinality(frame_b) > 0,
            "Neither connector frame_a nor frame_b of FixedTranslation object is connected");

          frame_b.r_0 = frame_a.r_0 +  Modelica.Mechanics.MultiBody.Frames.resolve1(frame_a.R, origin+u);

          if Connections.rooted(frame_a.R) then
            frame_b.R = Modelica.Mechanics.MultiBody.Frames.absoluteRotation(frame_a.R, R_theta);
            zeros(3) = frame_a.f + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_b.f);
            zeros(3) = frame_a.t + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_b.t) - cross(origin+u,
              frame_a.f);
          else
            frame_a.R = Modelica.Mechanics.MultiBody.Frames.absoluteRotation(frame_b.R, R_theta);
            zeros(3) = frame_b.f + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_a.f);
            zeros(3) = frame_b.t + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_a.t) + cross(
              Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, origin+u), frame_b.f);
          end if;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(
              extent={{-54,50},{54,-58}},
              lineColor={0,0,0},
              fillColor={170,213,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{36,36},{82,88}},
              color={0,0,0},
              thickness=0.5),
            Line(
              points={{-78,86},{-34,38}},
              color={0,0,0},
              thickness=0.5),
            Line(
              points={{-38,-44},{-86,-94}},
              color={0,0,0},
              thickness=0.5),
            Line(
              points={{32,-48},{94,-96}},
              color={0,0,0},
              thickness=0.5)}),                                          Diagram(
                coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                    100,100}})));
    end Node;
  end Components;

  class SID_File
    extends ExternalObject;

    function constructor
      input String fileName;
      output SID_File sid;
      external "C" sid = SIDFileConstructor_C(fileName)
      annotation(Include="#include \"ReadSID_C.h\"");

    end constructor;

    function destructor
      input SID_File sid;
      external "C" SIDFileDestructor_C(sid)
      annotation(Include="#include \"ReadSID_C.h\"");
    end destructor;

  end SID_File;

  package ExternalFunctions_C
    extends Modelica.Icons.FunctionsPackage;

     function getMass
      input EMBSlib.SID_File sid;
      output Real mass;
    external"C" mass = getMass(sid)
              annotation(Include="#include \"ReadSID_C.h\"");
     end getMass;

      function getM0
        input EMBSlib.SID_File sid;
        input String taylorName;
        input Integer nr;
        input Integer nc;
        output Real[nr,nc] m0;
        external "C" getM0(sid,taylorName,m0,nr,nc)
              annotation(Include="#include \"ReadSID_C.h\"");
      end getM0;

      function getM1
        input EMBSlib.SID_File sid;
        input String taylorName;
        input Integer nr;
        input Integer nq;
        input Integer nc;
        output Real[nr,nq,nc] m1;
        external "C" getM1(sid,taylorName,m1,nr,nq,nc)
              annotation(Include="#include \"ReadSID_C.h\"");
      end getM1;

      function getM0Node
        input EMBSlib.SID_File sid;
        input String taylorName;
        input Integer nodeIdx;
        input Integer nr;
        input Integer nc;
        output Real[nr,nc] m0;

        external "C" getM0Node(sid,taylorName,nodeIdx,m0,nr,nc)
              annotation(Include="#include \"ReadSID_C.h\"");
      end getM0Node;

      function getM1Node
        input EMBSlib.SID_File sid;
        input String taylorName;
        input Integer nodeIdx;
        input Integer nr;
        input Integer nq;
        input Integer nc;
        output Real[nr,nq,nc] m1;
        external "C" getM1Node(sid,taylorName,nodeIdx,m1,nr,nq,nc)
        annotation(Include="#include \"ReadSID_C.h\"");
      end getM1Node;

  end ExternalFunctions_C;

  package MatrixFunctions
      extends Modelica.Icons.FunctionsPackage;

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
            annotation(derivative=EMBSlib.MatrixFunctions.getTaylorFunction_d);
      end getTaylorFunction;

          function getGrMatrix
            input Integer nr;//(for Gr=3 for Ge=nq)
            input Integer nq;//nq
            input Integer nc;//3*nq
            input Real[nr,nc] Min;
            input Real[nq] q;
            output Real[nr, nr] Mout;
    protected
            Real[nr,nr] aux;
            Integer c=1;
          algorithm
            Mout := zeros(nr,nr);
            for i in 1:nq loop
              c :=(i - 1)*3;//column index
              aux[:,1] := Min[:,c+1];
              aux[:,2] := Min[:,c+2];
              aux[:,3] := Min[:,c+3];
              Mout := Mout+aux*q[i];
            end for;
          end getGrMatrix;

      function getTaylorFunction_d
            input Integer nr;
            input Integer nq;
            input Integer nc;
            input Real[nr,nc] M0;
            input Real[nr,nq,nc] M1;
            input Real[nq] q;
            input Real[nq] der_q;
            output Real[nr,nc] der_M;
    protected
            Real[nr,nc] aux;
      algorithm
            der_M := zeros(nr,nc);
            for i in 1:nq loop
              aux := M1[:,i,:]*q[i];
              der_M := der_M+aux;
            end for;
            der_M := M0 +der_M;

            annotation(derivative(order=2)=EMBSlib.MatrixFunctions.getTaylorFunction_2d);
      end getTaylorFunction_d;

      function getTaylorFunction_2d
            input Integer nr;
            input Integer nq;
            input Integer nc;
            input Real[nr,nc] M0;
            input Real[nr,nq,nc] M1;
            input Real[nq] q;
            input Real[nq] der_q;
            input Real[nq] der_2_q;
            output Real[nr,nc] der_2_M;
    protected
            Real[nr,nc] aux;
      algorithm
            der_2_M := zeros(nr,nc);
            for i in 1:nq loop
              aux := M1[:,i,:]*q[i];
              der_2_M := der_2_M+aux;
            end for;
            der_2_M := M0 +der_2_M;
      end getTaylorFunction_2d;

          function getGeMatrix
            input Integer nr;//(for Gr=3 for Ge=nq)
            input Integer nq;//nq
            input Integer nc;//3*nq
            input Real[nr,nc] Min;
            input Real[nq] q;
            output Real[nq, 3] Mout;
    protected
            Real[nq,3] aux;
            Integer c=1;
          algorithm
            Mout := zeros(nq,3);
            for i in 1:nq loop
              c :=(i - 1)*3;//column index
              aux[:,1] := Min[:,c+1];
              aux[:,2] := Min[:,c+2];
              aux[:,3] := Min[:,c+3];
              Mout := Mout+aux*q[i];
            end for;
          end getGeMatrix;
  end MatrixFunctions;

  model SimplePlate2 "Plate with force excitation"
    import FlexibleBodies;
    extends Modelica.Icons.Example;
    inner Modelica.Mechanics.MultiBody.World world(
      enableAnimation=true,
      n={0,0,-1},
      animateGravity=true,
      g=9.81,
      animateWorld=true)
      annotation (Placement(transformation(extent={{-112,10},{-92,30}},
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
      Real omega_tilde2[3,3] = omega_tilde*omega_tilde;


      // kinematic equations --> translation
      Real M_t[3] = mI*r_0_dd + transpose(mdCM_tilde)*phi_dd + transpose(Ct)*q_dd;
      Real k_omega_t[3] = 2*omega_tilde*transpose(Ct)*q_d + k_omega_t_2[:,1];
      Real k_omega_t_2[3,1] = omega_tilde*omega_tilde*mdCM_; // need the intermediate value to fix the dimensions
      Real f_comp[3] = M_t+k_omega_t;
      Real hd_t[3] = sum(identity(3)*ModalBody.nodes[i].f for i in 1:ModalBody.n_Nodes);
      Real residual_t[3] = f_comp - hd_t - ModalBody.frame_ref.f;  // this should be zero

      // kinematic equations --> rotation
      Real M_r[3] = mdCM_tilde*r_0_dd + J*phi_dd + transpose(Cr) *q_dd;
      Real k_omega_r[3] = Gr*phi_d + omega_tilde*J*phi_d;
      Real k_omega_r_1[3] = Gr*phi_d;
      Real k_omega_r_2[3] = omega_tilde*J*phi_d;
      Real t_comp[3] = M_r+k_omega_r;
      Real hd_r[3] = sum(ModalBody.nodes[i].r_0.*ModalBody.nodes[i].f+ModalBody.nodes[i].t for i in 1:ModalBody.n_Nodes);
      Real residual_r[3] = t_comp - hd_r - ModalBody.frame_ref.t;  // this should be zero

      //kinematic equations --> modal
      Real M_q[3] = Ct*r_0_dd + Cr*phi_dd + Me*q_dd;
      Real k_omega_q[3] = Ge*phi_d + Oe*OMega;
      Real k_omega_q_1[3] = Ge*phi_d;
      Real k_omega_q_2[3] = Oe*OMega;
      Real k_q[nq] = ksigma[:,1] + Ke*q +De*q_d;
      Real q_comp[nq] = M_q + k_omega_q + k_q;
      Real hde[nq] = phi*ModalBody.nodes[forceNode].f;

      //node position

      //disposition
      Real u[3]= ModalBody.nBody.u_out[forceNode, :];//forceNode
      //Real orig_M0[3,1] = ModalBody.nBody.modal.nodes[forceNode].origin;
      //Real orig_M1[3,nq,3] = ModalBody.nBody.modal.nodes[forceNode].origin;
      //Real orig [3,1] = EMBSlib.SID.MatrixFunctions.getTaylorFunction(3,nq,1,orig_M0,zeros(3,nq,1),q);
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


    Modelica.Mechanics.MultiBody.Joints.Revolute revFix(n(displayUnit="1") = {1,
        0,0}, w(fixed=true, start=0))
      annotation (Placement(transformation(extent={{-62,10},{-42,30}})));
  algorithm

  equation


    connect(world.frame_b, revFix.frame_a) annotation (Line(
        points={{-92,20},{-62,20}},
        color={95,95,95},
        thickness=0.5));
    connect(ModalBody.frame_ref, revFix.frame_b) annotation (Line(
        points={{-20,20},{-42,20}},
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
</html>"),   experiment(
        StopTime=10,
        Tolerance=1e-06,
        __Dymola_Algorithm="Dassl"),
      experimentSetupOutput,
      Commands(file="modelica://FlexibleBodies/Resources/Scripts/Settings.mos"
          "Settings",
        file="modelica://FlexibleBodies/Resources/Scripts/SimulateandPlotSimplePlate.mos"
          "Simulate and Plot Results",
        file="modelica://FlexibleBodies/Resources/Scripts/PlotSimplePlate.mos"
          "Plot Results"));
  end SimplePlate2;
  annotation (uses(Modelica(version="3.2.3")));
end EMBSlib;
