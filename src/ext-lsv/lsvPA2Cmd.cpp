#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <stdio.h>
#include "base/abc/abc.h"
#include "bdd/cudd/cudd.h"
#include "bdd/cudd/cuddInt.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"

static int Lsv_CommandSymBdd(Abc_Frame_t* pAbc, int argc, char** argv);
static int Lsv_CommandSymSat(Abc_Frame_t* pAbc, int argc, char** argv);
static int Lsv_CommandSymAll(Abc_Frame_t* pAbc, int argc, char** argv);

void init(Abc_Frame_t* pAbc)
{
  Cmd_CommandAdd(pAbc, "LSV", "lsv_sym_bdd", Lsv_CommandSymBdd, 1);
  Cmd_CommandAdd(pAbc, "LSV", "lsv_sym_sat", Lsv_CommandSymSat, 2);
  Cmd_CommandAdd(pAbc, "LSV", "lsv_sym_all", Lsv_CommandSymAll, 2);
}

void destroy(Abc_Frame_t* pAbc) {}

Abc_FrameInitializer_t frame_initializer = { init, destroy };

struct PackageRegistrationManager
{
  PackageRegistrationManager() { Abc_FrameAddInitializer(&frame_initializer); }
} lsvPackageRegistrationManager;

int Lsv_CommandSymBdd(Abc_Frame_t* pAbc, int argc, char** argv)
{
  Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
  Nm_Man_t* pManName = pNtk->pManName;
  int c;
  Extra_UtilGetoptReset();
  while ((c = Extra_UtilGetopt(argc, argv, "h")) != EOF)
  {
    switch (c)
    {
    case 'h':
      goto usage;
    default:
      goto usage;
    }
  }

  int outputIndex, input1index, input2index;
  if (argc == 3)
  {
    outputIndex = std::stoi(argv[1]);
    input1index = std::stoi(argv[2]);
    input2index = std::stoi(argv[3]);
  }
  else{
    goto usage;
  }

  if (!pNtk)
  {
    Abc_Print(-1, "Empty network.\n");
    return 1;
  }
  // Lsv_NtkPrintNodes(pNtk);
  if (!(Abc_NtkHasBdd(pNtk)))
  {
    printf("BDD is not found\n");
  }
  else
  {
    int i;
    bool symmetric = true;
    Abc_Obj_t* pPo = Abc_NtkPo(pNtk, outputIndex);
    Abc_Obj_t* pNode = Abc_ObjFanin0(pPo);
    Abc_Obj_t* pFanin = nullptr;
    DdManager* pDDman = (DdManager*)(pNtk->pManFunc);
    DdNode* poDDNode = (DdNode*)(pNode->pData);

    DdNode* pBDD1 = poDDNode;
    DdNode* pBDD2 = poDDNode;
    int numberOfPIs = Abc_NtkPiNum(pNtk);
    int numberOfPICombinations = 1 >> numberOfPIs;
    std::vector<std::vector<int>> prevBdd1patterns(1, std::vector<int>(0));
    std::vector<std::vector<int>> prevBdd2patterns(1, std::vector<int>(0));
    std::vector<std::vector<int>> bdd1patterns;
    std::vector<std::vector<int>> bdd2patterns;
    std::vector<DdNode*> bddNodesPrev1(1, pBDD1);
    std::vector<DdNode*> bddNodes1;
    std::vector<DdNode*> bddNodesPrev2(1, pBDD2);
    std::vector<DdNode*> bddNodes2;
    std::vector<int> assymetricPattern1;
    std::vector<int> assymetricPattern2;
    i = 0;
    Abc_ObjForEachFanin(pNode, pFanin, i)
    {
      DdNode* p11;
      DdNode* p12;
      DdNode* p21;
      DdNode* p22;
      for(int j = 0;j<bddNodesPrev1.size();++j){
        if(i == input1index){
          p11 = Cudd_Cofactor(pDDman, bddNodesPrev2[j], Cudd_Not(Cudd_bddIthVar(pDDman, i)));
          Cudd_Ref(p11);
          p22 = Cudd_Cofactor(pDDman, bddNodesPrev1[j], Cudd_bddIthVar(pDDman, i));
          Cudd_Ref(p22);
          bddNodes1.push_back(p11);
          bddNodes2.push_back(p22);
          std::vector<int> pattern11 = prevBdd1patterns[j];
          std::vector<int> pattern22 = prevBdd2patterns[j];
          pattern11.push_back(0);
          pattern22.push_back(1);
          bdd1patterns.push_back(pattern11);
          bdd2patterns.push_back(pattern22);

        }
        else if(i == input2index){
          p12 = Cudd_Cofactor(pDDman, bddNodesPrev2[j], Cudd_bddIthVar(pDDman, i));
          Cudd_Ref(p12);
          p21 = Cudd_Cofactor(pDDman, bddNodesPrev1[j], Cudd_Not(Cudd_bddIthVar(pDDman, i)));
          Cudd_Ref(p21);
          bddNodes1.push_back(p12);
          bddNodes2.push_back(p21);
          std::vector<int> pattern12 = prevBdd1patterns[j];
          std::vector<int> pattern21 = prevBdd2patterns[j];
          pattern12.push_back(1);
          pattern21.push_back(0);
          bdd1patterns.push_back(pattern12);
          bdd2patterns.push_back(pattern21);
        }
        else{
          p11 = Cudd_Cofactor(pDDman, bddNodesPrev1[j], Cudd_Not(Cudd_bddIthVar(pDDman, i)));
          Cudd_Ref(p11);
          p12 = Cudd_Cofactor(pDDman, bddNodesPrev1[j], Cudd_bddIthVar(pDDman, i));
          Cudd_Ref(p12);
          p21 = Cudd_Cofactor(pDDman, bddNodesPrev2[j], Cudd_Not(Cudd_bddIthVar(pDDman, i)));
          Cudd_Ref(p21);
          p22 = Cudd_Cofactor(pDDman, bddNodesPrev2[j], Cudd_bddIthVar(pDDman, i));
          Cudd_Ref(p22);
          bddNodes1.push_back(p11);
          bddNodes1.push_back(p12);
          bddNodes2.push_back(p21);
          bddNodes2.push_back(p22);
          std::vector<int> pattern11 = prevBdd1patterns[j];
          std::vector<int> pattern12 = prevBdd1patterns[j];
          std::vector<int> pattern21 = prevBdd2patterns[j];
          std::vector<int> pattern22 = prevBdd2patterns[j];
          pattern11.push_back(0);
          pattern12.push_back(1);
          pattern21.push_back(0);
          pattern22.push_back(1);
          bdd1patterns.push_back(pattern11);
          bdd1patterns.push_back(pattern12);
          bdd2patterns.push_back(pattern21);
          bdd2patterns.push_back(pattern22);
        }
      }
      for(int j = 0;j<bddNodesPrev1.size();++j){
        Cudd_RecursiveDeref(pDDman, bddNodesPrev1[j]);
        Cudd_RecursiveDeref(pDDman, bddNodesPrev2[j]);
      }

      if (i == (numberOfPIs - 1)){
        for (int j = 0;j < bddNodes1.size();++j){
          if(((bddNodes1[j] == Cudd_ReadOne(pDDman)) && (bddNodes2[j] == Cudd_ReadLogicZero(pDDman))) 
            || (((bddNodes2[j] == Cudd_ReadOne(pDDman)) && (bddNodes1[j] == Cudd_ReadLogicZero(pDDman))))) {
              symmetric = false;
              assymetricPattern1 = bdd1patterns[j];
              assymetricPattern2 = bdd2patterns[j];
              break;
          }
        }
      }
      bddNodesPrev1 = bddNodes1;
      bddNodesPrev2 = bddNodes2;
      bddNodes1.clear();
      bddNodes2.clear();
      prevBdd1patterns = bdd1patterns;
      prevBdd2patterns = bdd2patterns;
      bdd1patterns.clear();
      bdd2patterns.clear();
    }
    if(symmetric){
      std::cout << "symmetric\n";
    }
    else{
      std::string pattern1 = "";
      std::string pattern2 = "";
      for(int j = 0;j<numberOfPIs;++j){
        if(assymetricPattern1[j] == 0){
          pattern1 += '0';
        }
        else{
          pattern1 += '1';
        }
        if (assymetricPattern2[j] == 0) {
          pattern2 += '0';
        }
        else {
          pattern2 += '1';
        }
      }
      std::cout << "asymmetric\n";
      std::cout << pattern1 << "\n";
      std::cout << pattern2 << "\n";
    }
  }

usage:
  Abc_Print(-2, "usage: lsv_sym_bdd [-h]\n");
  Abc_Print(-2, "\t-h    : print the command usage\n");
  return 1;
}

int Lsv_CommandSymSat(Abc_Frame_t* pAbc, int argc, char** argv)
{
  return 0;
}

int Lsv_CommandSymAll(Abc_Frame_t* pAbc, int argc, char** argv)
{
  return 0;
}
