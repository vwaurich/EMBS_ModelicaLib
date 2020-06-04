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
      //Taylortest[numNodes] nodes;
      Taylortest J;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SID_DataStructure;

    record Taylortest
      Integer nrow;
      Integer ncol;
      Real[nrow,ncol] M0; //nrows,ncol
      //Real[nrow,nq,ncol] M1; //rows,q,cols
    end Taylortest;
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
      parameter Integer nrow;
      parameter Integer ncol;
      parameter Real[nrow,ncol] M0 = EMBSlib.SID.ParserFunctions.getTaylor(nrow);


      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Taylor;

    model SID_File

      parameter String fileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
      parameter EMBSlib.Types.SID_DataStructure sid = EMBSlib.SID.ParserFunctions.getSID_DataStructure(fileName) annotation(Evaluate=true);
      parameter EMBSlib.Types.Taylortest[sid.numNodes] nodes = EMBSlib.SID.ParserFunctions.getNodes(sid.numNodes)  annotation(Evaluate=true);
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
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

      function getSID_DataStructure
        input String fileName;
        output EMBSlib.Types.SID_DataStructure struc;
      protected
        Boolean endOfFile=false;
        String file;
        String line;
        Modelica.Utilities.Types.TokenValue token;
        Integer iline=1;
        Integer nextIndex,nextTokenIdx;
        Boolean found;
        EMBSlib.Types.Taylortest t;
      algorithm
         file :=Modelica.Utilities.Files.loadResource(fileName);

         // the first 2 integer tokens are numNodes and numModes
         (line, endOfFile) := Modelica.Utilities.Streams.readLine(file, 1);
         Modelica.Utilities.Streams.print(line);
         (struc.numNodes,found,nextTokenIdx) := EMBSlib.SID.ParserFunctions.getNextInteger(line,1);
         (struc.numModes,found,nextTokenIdx) := EMBSlib.SID.ParserFunctions.getNextInteger(line, nextTokenIdx);
         struc.J.nrow:=3;
         struc.J.ncol:=struc.numNodes;
         struc.J.M0 := EMBSlib.SID.ParserFunctions.getTaylor(struc.J.nrow,struc.J.ncol);

         //struc.nodes :=EMBSlib.SID.ParserFunctions.getNodes(3);

         while not found and not endOfFile loop
          (token, nextIndex) := Modelica.Utilities.Strings.scanToken(line);
          found := true;
          // read next line
          //(line, endOfFile) := Modelica.Utilities.Streams.readLine(fileName, iline);
         end while;
      end getSID_DataStructure;

      function getTaylor
        input Integer nrow;
        input Integer ncol;
        output Real[nrow,ncol] M0;
      algorithm
        for i in 1:nrow loop
          M0[i,1] :=i;
          M0[i,2] :=i+1;
        end for;
      end getTaylor;

      function testF
        output Integer i;
      algorithm
        i:=3;
      end testF;

      function getNodes
        input Integer size;
        output EMBSlib.Types.Taylortest[size] node;
      algorithm
        for i in 1:size loop
          node[i].nrow := 2;
          node[i].ncol := 2;
          node[i].M0 := EMBSlib.SID.ParserFunctions.getTaylor(2,2);
        end for;
      end getNodes;
    end ParserFunctions;
  end SID;
  annotation (uses(Modelica(version="3.2.3")));
end EMBSlib;
