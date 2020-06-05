within ;
package EMBSlib

  model EMBS_bodyExample
    Components.EMBS_Body eMBS_Body
      annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  equation
    connect(eMBS_Body.frame_a, world.frame_b) annotation (Line(
        points={{-60,10.2},{-60,10},{-80,10}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    annotation (Diagram(graphics));
  end EMBS_bodyExample;

  package Components
    model EMBS_Body
      parameter Integer numNodes = 9;

      Nodes nodes(numNodes=numNodes)
                  annotation (Placement(transformation(extent={{0,40},{20,60}})));

    Real              rShapes[numNodes,3];
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape[numNodes](r= rShapes,
        shapeType="sphere",
        length=0.1,
        width=0.1,
        height=0.1)
        annotation (Placement(transformation(extent={{80,60},{100,80}})));
      Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a
        annotation (Placement(transformation(extent={{-116,-14},{-84,18}})));
      Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_node[numNodes]
        annotation (Placement(transformation(extent={{84,-14},{116,18}})));
    algorithm
        for i in 1:numNodes loop
          rShapes[i,1] := nodes.origin[i].M0[1,1];
          rShapes[i,2] := nodes.origin[i].M0[2,1];
          rShapes[i,3] := nodes.origin[i].M0[3,1];
        end for;
    end EMBS_Body;

    model Nodes
      parameter Integer numNodes = 9;
      parameter String SIDfileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
       EMBSlib.Types.Taylor origin[numNodes];
       EMBSlib.Types.Taylor phi[numNodes];
       EMBSlib.Types.Taylor psi[numNodes];
       EMBSlib.Types.Taylor AP[numNodes];
       EMBSlib.Types.Taylor J;
       EMBSlib.SID_Data sid = EMBSlib.SID_Data(SIDfileName);

    algorithm
      J :=  EMBSlib.Functions.getTaylor(sid, "J");

      for i in 1:numNodes loop

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

      end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));
    end Nodes;

    model Node
      parameter Integer numNodes = 9;
      parameter String SIDfileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
      EMBSlib.Types.Taylor  J = EMBSlib.Functions.getTaylor(sid, "J");
      EMBSlib.SID_Data sid = EMBSlib.SID_Data(SIDfileName) annotation (Evaluate=true);

    algorithm
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));
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

    //remaining Taylor================================================
    //=========================================================

    pure function getTaylorOrder
      input EMBSlib.SID_Data sid;
      input String tName;
      output Integer order;
    external"C" order = getTaylorOrder(sid, tName) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getTaylorOrder;

    function getTaylorNrow
      input EMBSlib.SID_Data sid;
      input String tName;
      output Integer nrow;
    external"C" nrow = getTaylorNrow(sid, tName) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getTaylorNrow;

    function getTaylorNcol
      input EMBSlib.SID_Data sid;
      input String tName;
      output Integer ncol;
    external"C" ncol = getTaylorNcol(sid, tName) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getTaylorNcol;

    function getTaylorNq
      input EMBSlib.SID_Data sid;
      input String tName;
      output Integer nq;
    external"C" nq = getTaylorNq(sid, tName) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getTaylorNq;

    function getTaylorNqn
      input EMBSlib.SID_Data sid;
      input String tName;
      output Integer nqn;
    external"C" nqn = getTaylorNqn(sid, tName) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getTaylorNqn;

    function getTaylorStructure
      input EMBSlib.SID_Data sid;
      input String tName;
      output Integer nqn;
    external"C" nqn = getTaylorStructure(sid, tName) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getTaylorStructure;

    function getTaylorM0
      input EMBSlib.SID_Data sid;
      input String tName;
      input Integer r;
      input Integer c;
      input Integer dimR;
      input Integer dimC;
      output Real value;
      external "C" value = getTaylorM0(sid,tName,r,c, dimR, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getTaylorM0;

      function getTaylorM1
      input EMBSlib.SID_Data sid;
      input String tName;
      input Integer r;
      input Integer q;
      input Integer c;
      input Integer dimR;
      input Integer dimQ;
      input Integer dimC;
      output Real value;
      external "C" value = getTaylorM1(sid,tName,r,q,c, dimR,dimQ, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getTaylorM1;
  end ExternalFunctions;

  package Types
    record Taylor
      Integer order;
      Integer nrow;
      Integer ncol;
      Integer nq;
      Integer nqn;
      Integer structure;
      //Real[nrow,ncol] M0; //nrows,ncol
      //Real[nrow,nq,ncol] M1; //rows,q,cols
    end Taylor;

    record Matrix31
        Real M[3];
    end Matrix31;

    record SID_DataStructure
      Integer numNodes;
      Integer numModes;
      Modal modal;
      Taylor mdCM;
      Taylor J;
      Taylor Ct;
      Taylor Cr;
      Taylor Me;
      Taylor Gr;
      Taylor Ge;
      Taylor Oe;
      Taylor ksigma;
      Taylor Ke;
      Taylor De;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SID_DataStructure;

    record Taylortest
      Integer nrow;
      Integer ncol;
      //Real[nrow,ncol] M0; //nrows,ncol
      //Real[nrow,nq,ncol] M1; //rows,q,cols
    end Taylortest;

    record Modal
      Real mass;
      Integer nelastq;
      //Real[nelastq] freq;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Modal;

    record Nodes
    end Nodes;

    record Taylor2
      Integer nrow;
      Integer ncol;
      Real[nrow,ncol] M0; //nrows,ncol
      //Real[nrow,nq,ncol] M1; //rows,q,cols
    end Taylor2;
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

     function getTaylor
        input EMBSlib.SID_Data sid;
        input String tName;
        output EMBSlib.Types.Taylor taylor;
     algorithm
       taylor.order := EMBSlib.ExternalFunctions.getTaylorOrder(sid, tName)
                                                                           annotation(Evaluate=true);
       taylor.nrow := EMBSlib.ExternalFunctions.getTaylorNrow(sid, tName)
                                                                         annotation(Evaluate=true);
       taylor.ncol := EMBSlib.ExternalFunctions.getTaylorNcol(sid, tName)
                                                                         annotation(Evaluate=true);
       taylor.nq := EMBSlib.ExternalFunctions.getTaylorNq(sid, tName)
                                                                     annotation(Evaluate=true);
       taylor.nqn := EMBSlib.ExternalFunctions.getTaylorNqn(sid, tName)
                                                                       annotation(Evaluate=true);

       taylor.structure := EMBSlib.ExternalFunctions.getTaylorStructure(sid, tName)
                                                                                   annotation(Evaluate=true);
        for r0 in (1:taylor.nrow) loop
          for c0 in (1:taylor.ncol) loop
              //taylor.M0[r0,c0] := EMBSlib.ExternalFunctions.getTaylorM0(sid, tName,r0,c0,taylor.nrow,taylor.ncol);
      end for;
        end for;

      for r1 in (1:taylor.nrow) loop
        for q1 in (1:taylor.nq) loop
          for c1 in (1:taylor.ncol) loop
              //taylor.M1[r1,q1,c1] := EMBSlib.ExternalFunctions.getTaylorM1(sid, tName,r1,q1,c1,taylor.nrow,taylor.nq, taylor.ncol);
      end for;
      end for;
      end for;
     end getTaylor;

  end Functions;

  package SID
    model Taylor
      parameter String file;
      parameter Integer nodeId;
      parameter Integer nodeLineIdx = EMBSlib.SID.ParserFunctions.getLineIndexForNode(file,nodeId);
      parameter Integer nrow = EMBSlib.SID.ParserFunctions.getNodeNumRow(file,nodeLineIdx,nodeId);
      parameter Integer ncol = EMBSlib.SID.ParserFunctions.getNodeNumCol(file,nodeLineIdx,nodeId);
      parameter Real[nrow,ncol] origin_M0 = EMBSlib.SID.ParserFunctions.parseMatrix2(file,nodeLineIdx,nrow,ncol);

      Real[3] rShape;
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape(
        r=rShape,
        shapeType="sphere",
        length=0.1,
        width=0.1,
        height=0.1) annotation (Placement(transformation(extent={{-80,62},{-60,82}})));
      Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (Placement(transformation(extent={{84,-14},{116,18}})));
    algorithm
          rShape[1] := origin_M0[1,1];
          rShape[2] := origin_M0[2,1];
          rShape[3] := origin_M0[3,1];

          frame_a.r_0 := rShape;
          frame_a.R := Modelica.Mechanics.MultiBody.Frames.nullRotation();

    end Taylor;

    model SID_File

      parameter String fileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
      parameter String file = Modelica.Utilities.Files.loadResource(fileName);

      parameter EMBSlib.Types.SID_DataStructure sid = EMBSlib.SID.ParserFunctions.getSID_DataStructure(fileName) annotation(Evaluate=true);

      parameter Integer n = size(nodeLst,1);
      parameter Integer[:] nodeLst =  {1,2,12,26,78,84,107,143,149};

      Taylor nodes[n]( each file = file,nodeId=nodeLst) annotation (Placement(transformation(extent={{-80,40},{-60,60}})));


    //check out: https://github.com/modelica/ModelicaSpecification/issues/2282



    //funktoniert nicht:


      //parameter Integer[3] rows = {3,2,2};
      //parameter Integer[3] cols = {3,2,2};
      //parameter Real M[3,2,2] = {{{1}},{{1}},{{1}}};
      //parameter Real M1[3,3] = {{1,2,4},{3,4,5},{3,4,5}};
      //parameter Real M2[2,2] = {{11,13},{18,71}};
      //parameter Real M3[2,2] = {{21,24},{25,26}};

      //parameter EMBSlib.Types.Taylor2[3] nodes(nrow=rows, ncol=cols,M0 = {M1,M2,M3})



        annotation (Placement(transformation(extent={{-78,0},{-58,20}})),
                 Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end SID_File;

    package ParserFunctions

      function getNextInteger
        input String line;
        input Integer startIdx;
        output Integer i;
        output Boolean found = false;
        output Integer nextIndex = startIdx+1;
      protected
        Modelica.Utilities.Types.TokenValue token;
        Integer length = Modelica.Utilities.Strings.length(line);
      algorithm
         while
              (not found) loop
           (token, nextIndex) := Modelica.Utilities.Strings.scanToken(line,startIdx);
           if
             (token.tokenType == Modelica.Utilities.Types.TokenType.IntegerToken) then
             i := token.integer;
             found :=true;
           end if;
           if
             (nextIndex>=length) then
             i :=-1;
             break;
           end if;

         end while;

      end getNextInteger;

      function getRealValue "searches for the first real value after an equation mark. values can have the format x.xxxD+yy"
        input String line;
        input Integer startIdx;
        input String delimiter = "";
        output Boolean found = false;
        output Real value = 0.0;
      protected
        Modelica.Utilities.Types.TokenValue token;
        Integer length = Modelica.Utilities.Strings.length(line);
        Integer nextIdx = startIdx;
      algorithm
         Modelica.Utilities.Streams.print("LINE "+line);
         while
              (not found) loop
           (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
            //Modelica.Utilities.Streams.print("Token "+String(token.tokenType));
           if (token.tokenType == Modelica.Utilities.Types.TokenType.DelimiterToken) and token.string==delimiter then
             //Modelica.Utilities.Streams.print("Delimiter Token "+token.string);
             (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
             if (token.tokenType == Modelica.Utilities.Types.TokenType.RealToken) then
               //Modelica.Utilities.Streams.print("Got Real Token "+String(token.real));
               value := token.real;
               found := true;
               (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
               //Modelica.Utilities.Streams.print("Next Token "+String(token.tokenType)+" : "+token.string);
               if (token.tokenType == Modelica.Utilities.Types.TokenType.IdentifierToken) and token.string=="D" then
                 (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
                  //Modelica.Utilities.Streams.print("Next Next Token "+String(token.tokenType)+" : "+String(token.integer));
                 if
                   (token.tokenType == Modelica.Utilities.Types.TokenType.IntegerToken) then
                   value := value*10^token.integer;
                 end if;
               end if;
             end if;
           end if;
           if
             (nextIdx>=length) then
             break;
           end if;
         end while;

      end getRealValue;

      function getIntegerValue
        "searches for the first integer value after an equation mark"
        input String line;
        input Integer startIdx;
        input String delimiter = "=";
        output Boolean found = false;
        output Integer value = 0;
      protected
        Modelica.Utilities.Types.TokenValue token;
        Integer length = Modelica.Utilities.Strings.length(line);
        Integer nextIdx = startIdx;
      algorithm

         //Modelica.Utilities.Streams.print("getIntegerValue LINE "+line+" startIDx "+String(nextIdx));
         while
              (not found) loop
           (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
            //Modelica.Utilities.Streams.print("Token "+String(token.tokenType)+"  "+token.string);
           if (token.tokenType == Modelica.Utilities.Types.TokenType.DelimiterToken) and token.string==delimiter then
             //Modelica.Utilities.Streams.print("Delimiter Token "+token.string);
             (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
             if (token.tokenType == Modelica.Utilities.Types.TokenType.IntegerToken) then
               //Modelica.Utilities.Streams.print("Got Int Token "+String(token.integer));
               value := token.integer;
               found := true;
             end if;
           end if;
           if
             (nextIdx>=length) then
             break;
           end if;
         end while;
      end getIntegerValue;

      function getSID_DataStructure
        input String file;
        output EMBSlib.Types.SID_DataStructure struc;
      protected
        Boolean endOfFile=false;
        Boolean found = false;
        String line;
        Modelica.Utilities.Types.TokenValue token;
        Integer lineIdx=1;
        Integer nextIndex,nextTokenIdx;
        Boolean foundPart=false;
      algorithm
         //file :=Modelica.Utilities.Files.loadResource(fileName);

         // the first 2 integer tokens are numNodes and numModes
         (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
         (struc.numNodes,found,nextTokenIdx) := EMBSlib.SID.ParserFunctions.getNextInteger(line,1);
         (struc.numModes,found,nextTokenIdx) := EMBSlib.SID.ParserFunctions.getNextInteger(line, nextTokenIdx);
         //struc.nodes := EMBSlib.SID.ParserFunctions.getNodes(2);

         //get the modal information
         //(struc.modal.mass,struc.modal.nelastq,lineIdx) := parseModalData(file,lineIdx) annotation(Evaluate=true);
         //struc.modal.freq := parseModalFrequencies(file,lineIdx,struc.modal.nelastq);

         //traverse the node data

         //traverse the additional data
         //(struc.mdCM,lineIdx) := parseTaylorData( file,lineIdx, "mdCM");
         //(struc.J,lineIdx) := parseTaylorData( file,lineIdx, "J");
         //(struc.Ct,lineIdx) := parseTaylorData( file,lineIdx, "Ct");
         //(struc.Cr,lineIdx) := parseTaylorData( file,lineIdx, "Cr");
         //(struc.Me,lineIdx) := parseTaylorData( file,lineIdx, "Me");
         //(struc.Gr,lineIdx) := parseTaylorData( file,lineIdx, "Gr");
         //(struc.Ge,lineIdx) := parseTaylorData( file,lineIdx, "Ge");
         //(struc.Oe,lineIdx) := parseTaylorData( file,lineIdx, "Oe");
         //(struc.ksigma,lineIdx) := parseTaylorData( file,lineIdx, "ksigma");
         //(struc.Ke,lineIdx) := parseTaylorData( file,lineIdx, "Ke");
         //(struc.De,lineIdx) := parseTaylorData( file,lineIdx, "De");
      end getSID_DataStructure;

      function getSID_DataStructureNodes
        input String fileName;
        input Integer size;
        input EMBSlib.Types.SID_DataStructure struc;
        output EMBSlib.Types.Taylor[size] tOut;
      protected
        Boolean endOfFile=false;
        Boolean found = false;
        String file;
        String line;
        Modelica.Utilities.Types.TokenValue token;
        Integer lineIdx=1;
        Integer nextIndex,nextTokenIdx;
        Boolean foundPart=false;
        EMBSlib.Types.Taylor[size] t;
      algorithm
         file :=Modelica.Utilities.Files.loadResource(fileName);

         for i in 1:size loop
              //tOut[i].order := 1;
              //tOut[i].nrow := 2;
              //tOut[i].ncol := 3;
              //tOut[i].nq := 4;
              //tOut[i].nqn := 5;
              //tOut[i].structure := 6;
              //tOut[i].M0 := parseMatrix2(file, 1,tOut[i].nrow,tOut[i].ncol);
         end for;



      end getSID_DataStructureNodes;

      function getTaylor
        input Integer nrow;
        input Integer ncol;
        output Real[nrow,ncol] M0;
      algorithm
        for i in 1:nrow loop
          M0[i,1] :=i;
          M0[i,1] :=i+1;
        end for;
      end getTaylor;

      function getNodes
        input Integer size;
        output EMBSlib.Types.Taylor[size] node;
      protected
        Integer nrow;
        Integer ncol;
      algorithm
        for i in 1:size loop
          node[i].order := 1;
          node[i].nrow := 2+i;
          node[i].ncol := 2+i;
          node[i].nq := i;
          node[i].nqn := 3;
          node[i].structure := 2;
         // node[i].M0 := EMBSlib.SID.ParserFunctions.getTaylor(2,2);
        end for;
      end getNodes;

      function findObject
        input String file;
        input Integer lineIdxIn;
        input String objectName;
        output Boolean found = false;
        output Integer lineIdxOut;
      protected
        Boolean endOfFile = false;
        Boolean foundModal = false;
        Integer nextTokenIdx = 1;
        Integer lineLength;
        Integer lineIdx = max(1,lineIdxIn);
        String line;
        Modelica.Utilities.Types.TokenValue token;
      algorithm
        while not found and not endOfFile loop
          (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
          (found, nextTokenIdx) := EMBSlib.SID.ParserFunctions.checkForIdentifier(line,objectName);
          if (not found) then
            lineIdx := lineIdx+1;
          end if;
        end while;

        lineIdxOut := lineIdx;
      end findObject;

      function findObject2
        input String file;
        input Integer lineIdxIn;
        input String objectName1;
        input String objectName2;
        output Boolean found = false;
        output Integer lineIdxOut;
      protected
        Boolean endOfFile = false;
        Boolean foundModal = false;
        Integer nextTokenIdx = 1;
        Integer lineLength;
        Integer lineIdx = lineIdxIn;
        String line;
        Modelica.Utilities.Types.TokenValue token;
      algorithm
        while not found and not endOfFile loop
          (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
          (found, nextTokenIdx) := EMBSlib.SID.ParserFunctions.checkFor2Identifier(line,objectName1, objectName2);
          if (not found) then
            lineIdx := lineIdx+1;
          end if;
        end while;

        lineIdxOut := lineIdx;
      end findObject2;

      function findModal
        input String file;
        input Integer lineIdxIn;
        output Boolean foundModal = false;
        output Integer lineIdxOut;
      protected
        Boolean endOfFile = false;
        Integer nextIdx = 1;
        Integer lineLength;
        Integer lineIdx = lineIdxIn;
        String line;
        Modelica.Utilities.Types.TokenValue token;
      algorithm
        while not foundModal and not endOfFile loop
          (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
          (foundModal, nextTokenIdx) := EMBSlib.SID.ParserFunctions.checkFor2Identifier(line,"new","modal");
          lineIdx := lineIdx+1;
          if (foundModal) then
            Modelica.Utilities.Streams.print("found modal in "+line);
          end if;
        end while;

        lineIdxOut := lineIdx;
      end findModal;

      function checkForMatrix2Entry
        input String line;
        input String matrixIdent "for example \"m\"";
        output Boolean found = false;
        output Integer idx1;
        output Integer idx2;
        output Real value;
      protected
        Integer nextIdx = 1;
        Integer lineLength;
        Modelica.Utilities.Types.TokenValue token;
      algorithm
          lineLength := Modelica.Utilities.Strings.length(line);
          found := false;
          while not found and (nextIdx<lineLength) loop
            (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
             Modelica.Utilities.Streams.print("Token: "+String(token.tokenType)+"  "+token.string +" in line "+line);

            if (token.tokenType == Modelica.Utilities.Types.TokenType.IdentifierToken) then
              if (token.string == matrixIdent) then
                Modelica.Utilities.Streams.print("found matrix ident "+matrixIdent);
                (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
                found := true;
                Modelica.Utilities.Streams.print("found proper index ");
                nextIdx := nextIdx-1;
                (found,idx1) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,nextIdx,"(");
                (found,idx2) :=EMBSlib.SID.ParserFunctions.getIntegerValue(line,nextIdx, ",");
                (found,value) :=EMBSlib.SID.ParserFunctions.getRealValue(line, nextIdx,"=");
              else
                found := false;
              end if;
            end if;
          end while;
      end checkForMatrix2Entry;

      function checkFor2Identifier
        input String line;
        input String id1;
        input String id2;
        output Boolean found = false;
        output Integer nextTokenId;
      protected
        Integer nextIdx = 1;
        Integer lineLength;
        Modelica.Utilities.Types.TokenValue token;

      algorithm
          lineLength := Modelica.Utilities.Strings.length(line);
          while not found and nextIdx<lineLength loop
            (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
            if (token.tokenType == Modelica.Utilities.Types.TokenType.IdentifierToken) then
              if (token.string == id1) then
                Modelica.Utilities.Streams.print("found first token "+id1);
                (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
                if (token.tokenType == Modelica.Utilities.Types.TokenType.IdentifierToken) then
                  if (token.string == id2) then
                    Modelica.Utilities.Streams.print("found second token "+id2);
                    found :=true;
                  end if;
                end if;
              end if;
            end if;
          end while;
      end checkFor2Identifier;

      function checkForIdentifier
        input String line;
        input String id1;
        output Boolean found = false;
        output Integer nextTokenId;
      protected
        Integer nextIdx = 1;
        Integer lineLength;
        Modelica.Utilities.Types.TokenValue token;

      algorithm
          lineLength := Modelica.Utilities.Strings.length(line);
          while not found and nextIdx<lineLength loop
            (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
            if (token.tokenType == Modelica.Utilities.Types.TokenType.IdentifierToken) then
              if (token.string == id1) then
                (token, nextIdx) := Modelica.Utilities.Strings.scanToken(line,nextIdx);
                found :=true;
              end if;
            end if;
          end while;
      end checkForIdentifier;

      function parseModalData
        input String file;
        input Integer lineIdxIn;
        output Real mass = 1.4;
        output Integer nelastq = 2;
        output Integer lineIdxOut;
      protected
        Integer lineIdx;
        String line;
        Boolean endOfFile;

        Boolean found;
      algorithm
        (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdxIn,"refmod");
        if
          (found) then

            (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdxIn,"mass");
            if (found) then
              (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
              (found,mass) := EMBSlib.SID.ParserFunctions.getRealValue(line,1,"=");
            end if;

            (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdxIn,"nelastq");
            if (found) then
              (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
              (found,nelastq) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1);
            end if;
        end if;
        lineIdxOut := lineIdx;
      end parseModalData;

      function parseModalFrequencies
        input String file;
        input Integer lineIdxIn;
        input Integer numFreqs;
        output Real[numFreqs] freqs;
        output Integer lineIdxOut;
      protected
        Integer lineIdx;
        String line;
        Boolean endOfFile;
        Boolean found;
      algorithm
        (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdxIn,"ielastq");
        for i in 1: numFreqs loop
          (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
          (found,freqs[i]) := EMBSlib.SID.ParserFunctions.getRealValue(line,1,":");
          lineIdx := lineIdx+1;
        end for;

        lineIdxOut := lineIdx;
      end parseModalFrequencies;

      function parseTaylorData
        input String file;
        input Integer lineIdxIn;
        input String taylorName;
        output EMBSlib.Types.Taylor t;
        output Integer lineIdxOut;
      protected
        Integer lineIdx;
        String line;
        Boolean endOfFile;
        Boolean found;

        Integer idx1,idx2;
        Real value;
      algorithm
         //Modelica.Utilities.Streams.print("Search taylor in "+String(lineIdxIn));
        (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdxIn,taylorName);
        (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);

         //Modelica.Utilities.Streams.print(String(found)+"Found taylor "+taylorName+" in "+line+String(lineIdx));

        if (found) then

          (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"order");
          if
            (found) then
            (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
            (found,t.order) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1);
              lineIdx := lineIdx+1;
          end if;

          (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"nrow");
          if
            (found) then
            (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
            (found,t.nrow) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1);
              lineIdx := lineIdx+1;
          end if;

          (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"ncol");
          if
            (found) then
            (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
            (found,t.ncol) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1);
              lineIdx := lineIdx+1;
          end if;

          (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"nq");
          if
            (found) then
            (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
            (found,t.nq) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1);
            lineIdx := lineIdx+1;
          end if;

          (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"nqn");
          if
            (found) then
            (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
            (found,t.nqn) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1);
            lineIdx := lineIdx+1;
          end if;

          (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"structure");
          if
            (found) then
            (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
            (found,t.structure) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1);
            lineIdx := lineIdx+1;
          end if;

        end if;
        //Modelica.Utilities.Streams.print("GOT Taylor"+taylorName+" order: "+String(t.order));


        //parse matrix entries
        //found :=true;
        //while found and not endOfFile loop
        //  (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"m0");
        //  (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
        //  (found, idx1,idx2,value) := EMBSlib.SID.ParserFunctions.checkForMatrix2Entry(line,"m0");
        //  //Modelica.Utilities.Streams.print(String(found)+" Found Matrix entry m0: "+String(idx1)+" "+String(idx2)+"  "+String(value));
        //  lineIdx :=lineIdx + 1;
        //end while;
        //lineIdxOut  := lineIdx;

      end parseTaylorData;

      function parseMatrix2
        input String file;
        input Integer lineIdxIn;
        input Integer nrow;
        input Integer ncol;
        output Real [nrow,ncol] M0;
      protected
        Integer lineIdx=lineIdxIn;
        Integer lineIdxEnd;
        String line;
        Boolean endOfFile;
        Boolean found;
        Integer idx1,idx2;
        Real value;
      algorithm
        found :=true;
        Modelica.Utilities.Streams.print("parseMAtrix2 "+String(nrow)+": "+String(ncol)+" at line "+String(lineIdx));
        while found and not endOfFile loop
          (found,lineIdxEnd) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"end");
          (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"m0");
            Modelica.Utilities.Streams.print("next end  "+String(lineIdxEnd)+" and next m0 "+String(lineIdx));


          if
            (lineIdxEnd<lineIdx) then
            endOfFile :=true;
          else
           (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
           (found, idx1,idx2,value) := EMBSlib.SID.ParserFunctions.checkForMatrix2Entry(line,"m0");
            Modelica.Utilities.Streams.print(String(found)+" Found Matrix entry m0: "+String(idx1)+" "+String(idx2)+"  "+String(value));
            M0[idx1,idx2] := value;
          end if;
            lineIdx :=lineIdx + 1;
        end while;
      end parseMatrix2;

      function getNodeNumRow
        input String file;
        input Integer lineIdxIn;
        input Integer nodeId;
        output Integer numRow;
      protected
        Boolean found=false;
        Boolean endOfFile;
        Integer lineIdx=lineIdxIn;
        Integer nextTokenIdx;
        Integer idx;
        String line;
      algorithm
        while not found and not endOfFile loop
          (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
          (found, nextTokenIdx) := EMBSlib.SID.ParserFunctions.checkFor2Identifier(line,"new","node");
          if (found) then
            (found,idx) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1,"=");
            found := idx==nodeId;
            //correct node identifier found
            if
              (found) then
              (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"nrow");
               if
                 (found) then
                 (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
                 (found,numRow) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1, "=");
                end if;
            end if;
          end if;
          lineIdx := lineIdx+1;
        end while;
      end getNodeNumRow;

      function getNodeNumCol
        input String file;
        input Integer lineIdxIn;
        input Integer nodeId;
        output Integer numRow;
      protected
        Boolean found=false;
        Boolean endOfFile;
        Integer lineIdx=lineIdxIn;
        Integer nextTokenIdx;
        Integer idx;
        String line;
      algorithm
        while not found and not endOfFile loop
          (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
          (found, nextTokenIdx) := EMBSlib.SID.ParserFunctions.checkFor2Identifier(line,"new","node");
          if (found) then
            (found,idx) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1,"=");
            found := idx==nodeId;
            //correct node identifier found
            if
              (found) then
              (found,lineIdx) := EMBSlib.SID.ParserFunctions.findObject(file,lineIdx,"ncol");
               if
                 (found) then
                 (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
                 (found,numRow) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1, "=");
                end if;
            end if;
          end if;
          lineIdx := lineIdx+1;
        end while;
      end getNodeNumCol;

      function getLineIndexForNode
        input String file;
        input Integer nodeId;
        output Integer lineIdxOut;
      protected
        Boolean found=false;
        Boolean endOfFile;
        Integer lineIdx=1;
        Integer nextTokenIdx;
        Integer idx;
        String line;
      algorithm
        while not found and not endOfFile loop
          (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, lineIdx);
          (found, nextTokenIdx) := EMBSlib.SID.ParserFunctions.checkFor2Identifier(line,"new","node");
          if (found) then
            (found,idx) := EMBSlib.SID.ParserFunctions.getIntegerValue(line,1,"=");
            found := idx==nodeId;
            if (found) then
              lineIdxOut := lineIdx;
            end if;
          end if;
          lineIdx := lineIdx+1;
        end while;
      end getLineIndexForNode;
    end ParserFunctions;
  end SID;
  annotation (uses(Modelica(version="3.2.3")));
end EMBSlib;
