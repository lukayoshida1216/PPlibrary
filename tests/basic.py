from penpy.base import *
import unittest

class ConfigTestCase(unittest.TestCase):

    def test_config_exists(self): #それぞれの関数,メソッドのテスト
        self.assertIsNotNone(default_config)

    def test_config_content_type(self):
        for rtype in default_config:
            self.assertTrue("reactions" in default_config[rtype])
            self.assertTrue("species" in default_config[rtype])
            self.assertTrue(len(default_config[rtype]["reactions"]) > 0)

class ReactionFactoryTestCase(unittest.TestCase):
    def setUp(self):
        self.rfact = ReactionFactory_temp(default_config) 
    
    def test_Exo_getReactionList(self):
        N = Species("N",5,1.0, 10.0)
        edge3 = Edge("Exo",None,[N],[]) 
        rlist = self.rfact.getReactionList(edge3)
        self.assertTrue(len(rlist[1])==1)
        self.assertTrue(rlist[1][0].reactants[0] is N)  #反応物がNである
        self.assertTrue(rlist[1][0].count[N] == -1)   
        
    def test_PEN_getReactionList(self):  #例外処理（Inhibitor,self-start?)
        G1 = Species("G1",5,1.0, 10.0)    #名前、拡散係数、安定性
        A = Species("A",5,1.0, 10.0) 
        edge1 = Edge("PEN",G1, [A], [A]) 
        rlist = self.rfact.getReactionList(edge1)
        self.assertTrue(len(rlist[1])==11) #reactionの長さ
        self.assertTrue(rlist[1][0].reactants[0] is A) #TO DO
        self.assertTrue(rlist[1][0].reactants[1] is G1)
        self.assertTrue(rlist[1][0].products[0] is G1_in)
        self.assertTrue(rlist[1][1].reactants[0] is A)
        self.assertTrue(rlist[1][1].reactants[1] is G1)
        self.assertTrue(rlist[1][1].products[0] is G1_out)
        self.assertTrue(rlist[1][2].reactants[0] is A)
        self.assertTrue(rlist[1][2].reactants[1] is G1_out)
        self.assertTrue(rlist[1][2].products[0] is G1_both)
        self.assertTrue(rlist[1][3].reactants[0] is A)
        self.assertTrue(rlist[1][3].reactants[1] is G1_in)
        self.assertTrue(rlist[1][3].products[0] is G1_both)
        self.assertTrue(rlist[1][4].reactants[0] is G1_in)
        self.assertTrue(rlist[1][4].products[0] is A)
        self.assertTrue(rlist[1][4].products[1] is G1)
        self.assertTrue(rlist[1][5].reactants[0] is G1_out)
        self.assertTrue(rlist[1][5].products[0] is A)
        self.assertTrue(rlist[1][5].products[1] is G1)
        self.assertTrue(rlist[1][6].reactants[0] is G1_both)
        self.assertTrue(rlist[1][6].products[0] is A)
        self.assertTrue(rlist[1][6].products[1] is G1_out)
        self.assertTrue(rlist[1][7].reactants[0] is G1_both)
        self.assertTrue(rlist[1][7].products[0] is A)
        self.assertTrue(rlist[1][7].products[1] is G1_in)
        self.assertTrue(rlist[1][8].reactants[0] is G1_in)
        self.assertTrue(rlist[1][8].products[0] is G1_ext)
        self.assertTrue(rlist[1][9].reactants[0] is G1_both)
        self.assertTrue(rlist[1][9].products[0] is G1_ext)
        self.assertTrue(rlist[1][9].products[1] is A)
        self.assertTrue(rlist[1][10].reactants[0] is G1_ext)
        self.assetTrue(rlist[1][10].products[0] is G1_both)
        #self.assertTrue
        self.assertTrue(rlist[1][0].count[N] == -1) #TO DO
   

if __name__ == "__main__":  
    unittest.main() #テストスクリプトのコマンドライン用インターフェースを提供


    #PEN test
    rfact = ReactionFactory_temp(base_config)
    G1 = Species("G1",5,1.0, 10.0)    #名前、拡散係数、安定性  ❓ライブラリの中のclassの使い方
    A = Species("A",5,1.0, 10.0) 
    edge0 = Edge("PEN",G1, [A], [A]) #
    
    print(rfact.getReactionList(edge0))    #できた？
