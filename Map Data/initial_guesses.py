
import numpy as np

def spectrum_interval(x,y, lb, ub):
    lowerx = np.where(x>lb)
    upperx = np.where(x<ub)
    interx = np.intersect1d(lowerx,upperx)
    x = x[interx]  
    y = y[interx]
    return(x,y)

def main(material_type, xs, ys, second_pass, sample):
    if material_type == 'BP':
       
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)  

        I01 = max(ys)/25 # magnitude
        Gamma1 = 3 # FWHM
        cen1 = 362 # peak centre
        
        I02 = max(ys)/6 # magnitude
        Gamma2 = 3 # FWHM
        cen2 = 439
    
        I03 = max(ys) # magnitude
        Gamma3 = 3 # FWHM
        cen3 = 467
    

        bnds = ((0, I03*5), (cen1-15,cen1+15), (0,np.inf),
                (0, I03*5), (cen2-15,cen2+15), (0,np.inf),
                (0, I03*5), (cen3-15,cen3+15), (0,np.inf),
            (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3,  m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]  
            
    elif material_type == 'diamond':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 4 # FWHM
        cen1 = 1332
        
        # strict bound for singal lor, 
        # intensity should not be larger than maximum signal
        bnds = ((0, max(ys)), (cen1-10,cen1+10), (1.1,20), 
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]        


    elif material_type == 'LiP8':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 30 # FWHM
        cen1 = 286
        
        I02 = max(ys) # magnitude
        Gamma2 = 49 # FWHM
        cen2 = 344 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 7 # FWHM
        cen3 = 360  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 57 # FWHM
        cen4 = 394 #1560-1630
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 10 # FWHM
        cen5 = 437 #1560-1630
        
        I06 = max(ys)/2 # magnitude
        Gamma6 = 8 # FWHM
        cen6 = 464 #1560-1630
        
        bnds = ((0, max(ys)*2), (cen1-10,cen1+10), (2,100),
                (0, max(ys)*2), (cen2-10,cen2+10), (2,100),
                (0, max(ys)*2), (cen3-10,cen3+10), (2,100),
                (0, max(ys)*2), (cen4-10,cen4+10), (2,100),
                (0, max(ys)*2), (cen5-10,cen5+10), (2,100),
                (0, max(ys)*2), (cen6-10,cen6+10), (2,100),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
               'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'PNR':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 4 # FWHM
        cen1 = 362
        
        I02 = max(ys) # magnitude
        Gamma2 = 6 # FWHM
        cen2 = 439 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 5 # FWHM
        cen3 = 468  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 5 # FWHM
        cen4 = 521 #1560-1630
        
        
        bnds = ((0, max(ys)*2), (cen1-10,cen1+10), (2,100),
                (0, max(ys)*2), (cen2-10,cen2+10), (2,100),
                (0, max(ys)*2), (cen3-10,cen3+10), (2,100),
                (0, max(ys)*2), (cen4-10,cen4+10), (2,100),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'LFP':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        # LFP Peak
        I01 = max(ys)/5 # magnitude
        Gamma1 = 10 # FWHM
        cen1 = 950
        
        # D' Peak
        I02 = max(ys) # magnitude
        Gamma2 = 50 # FWHM
        cen2 = 1160 

        # D Peak
        I03 = max(ys) # magnitude
        Gamma3 = 50 # FWHM
        cen3 = 1310  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        # A Peak?
        I04 = max(ys)/5 # magnitude
        Gamma4 = 10 # FWHM
        cen4 = 1500 #1560-1630
        
        # G Peak        
        I05 = max(ys)/5 # magnitude
        Gamma5 = 10 # FWHM
        cen5 = 1600 #1560-1630
        
        
        bnds = ((0, max(ys)*2), (cen1-15,cen1+15), (1,50),
                (0, max(ys)*2), (cen2-30,cen2+30), (1,200),
                (0, max(ys)*2), (cen3-20,cen3+20), (1,200),
                (0, max(ys)*2), (cen4-50,cen4+50), (1,200),
                (0, max(ys)*2), (cen5-20,cen5+20), (1,200),
                (m-0.1,m+0.1), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'BG Graphite':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        I01 = max(ys) # magnitude
        Gamma1 = 20 # FWHM
        cen1 = 1315#1500
        
        I02 = max(ys) #y1590 # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1605
        
        I03 = max(ys) # magnitude
        Gamma3 = 17 # FWHM
        cen3 = 1576#1350
        q = -4
        
        bnds = ((0.1, I01*5), (cen1-20,cen1+20), (0,250),
                (0.1, I02*5), (cen2-10,cen2+10), (0,180), 
                (0.1, I03*5), (cen3-10,cen3+10), (0,30),
                (-500,-0.1), (m-0.01,m+0.01), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'q', 'm', 'c']  
            
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i] 

    elif material_type == 'SC Graphite':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        I01 = max(ys) # magnitude
        Gamma1 = 20 # FWHM
        cen1 = 1330#1500
        
        I02 = max(ys) #y1590 # magnitude
        Gamma2 = 17 # FWHM
        cen2 = 1585
        q = -4
        
        I03 = max(ys) # magnitude
        Gamma3 = 17 # FWHM
        cen3 = 2650#1350
        
        bnds = ((0.1, I01*5), (cen1-30,cen1+30), (0,250),
                (0.1, I02*5), (cen2-20,cen2+20), (0,180), 
                (0.1, I03*5), (cen3-60,cen3+60), (0,80),
                (-500,-0.1), (m-0.01,m+0.01), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'q', 'm', 'c']  
            
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i] 
            
    elif material_type == 'graphite full':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        I01 = max(ys)/10 # magnitude
        Gamma1 = 20 # FWHM
        cen1 = 1350#1500
        
        I02 = max(ys)/10 #y1590 # magnitude
        Gamma2 = 30 # FWHM
        cen2 = 1575
        
        I03 = max(ys) # magnitude
        Gamma3 = 3 # FWHM
        cen3 = 2440#1350
        
        I04 = max(ys)/10 #y1590 # magnitude
        Gamma4 = 30 # FWHM
        cen4 = 2685#1450
        
        I05 = max(ys)/10# y1590 # magnitude
        Gamma5 = 30# FWHM
        cen5 = 2720
        
        bnds = ((0.1, I01*5), (cen1-20,cen1+20), (0,np.inf),
                (0.1, I02*5), (cen2-20,cen2+20), (0,np.inf), 
                (0.1, I03*5), (cen3-20,cen3+20), (0,250),
                (0.1, I04*5), (cen4-20,cen4+20), (0,250),
                (0.1, I05*5), (cen5-20,cen5+20), (0,np.inf), #15
                (m-0.01,m+0.01), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, I04, cen4, Gamma4,
                I05, cen5, Gamma5, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'm', 'c']  
            
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i] 
            
    elif material_type == 'Graphite KC24':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        q = -4
        
        I01 = max(ys)/3
        Gamma1 = 40 # FWHM
        cen1 = 1585#1500
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1605
        
        bnds = ((0, I01*5), (cen1-10,cen1+10), (0,np.inf),
                (0, I02*5), (cen2-10,cen2+10), (0,np.inf), #15
            (-50,0), (m-0.001,m+0.001), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'q', 'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]  
            
    elif material_type == 'PtG20-5c':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        q = -4
        
        I01 = max(ys)/5
        Gamma1 = 40 # FWHM
        cen1 = 1350#1500
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1580
        
        I03 = max(ys)/10 # magnitude
        Gamma3 = 160 # FWHM
        cen3 = 1620
        
        I04 = max(ys)/2
        Gamma4 = 50
        cen4 = 2700
        
        bnds = ((0, I01*2), (cen1-15,cen1+15), (0,np.inf),
                (0, I02*5), (cen2-15,cen2+15), (0,np.inf),
                (0, I03*5), (cen3-7,cen3+7), (0,np.inf),
                (0, I04*5), (cen4-20,cen4+20), (0,np.inf),
            (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, I04, cen4, Gamma4,  m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i] 
            
    elif material_type == 'Ion Etched':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        y1590 = ys[np.where(xs>1590)[0][-1]] # finds the y value closest to
                                             # x = 1590
        
        I01 = y1590 #max(ys) # magnitude
        Gamma1 = 160 # FWHM
        cen1 = 1340#1500
        
        I02 = y1590 # magnitude
        Gamma2 = 30 # FWHM
        cen2 = 1550
        
        I03 = y1590 # magnitude
        Gamma3 = 160 # FWHM
        cen3 = 1600
        
        bnds = ((0, I01*2), (cen1-15,cen1+15), (0,np.inf),
                (0, I02*5), (cen2-15,cen2+15), (0,np.inf),
                (0, I03*5), (cen3-15,cen3+15), (0,np.inf),
            (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3,  m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]   
            
            
            
    elif material_type == 'ie2':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        q = -40
        
        y1590 = ys[np.where(xs>1590)[0][-1]] # finds the y value closest to
                                             # x = 1590
        
        I01 = y1590 #max(ys) # magnitude
        Gamma1 = 160 # FWHM
        cen1 = 1340#1500
        
        I02 = y1590 # magnitude
        Gamma2 = 160 # FWHM
        cen2 = 1600
        
        bnds = ((0, I01*2), (cen1-55,cen1+55), (0,np.inf),
                (0, I02*5), (cen2-55,cen2+55), (0,np.inf), #15
            (-50,0), (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'q', 'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]  
        
    elif material_type == 'Lasered Diamond Undoped':     
        
        y2 = np.mean(ys[:20])
        y1 = np.mean(ys[-20:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        y1580 = ys[np.where(xs>1580)[0][-1]]
        #print(y1580)
        
        # I01 = max(ys) # magnitude
        # Gamma1 = 3 # FWHM
        # cen1 = 1332#1500
        
        # I02 = 1
        # Gamma2 = 2 # FWHM
        # cen2 = 1350
        # ming2 = 1
        
        # I03 = 1
        # Gamma3 = 2# FWHM
        # cen3 = 1590
        # q = -6
        
        #if second_pass == True:
            
        I01 = max(ys)/10 # magnitude
        Gamma1 = 2 # FWHM
        cen1 = 1332#1500
        
        I02 = y1580*1.1#max(ys)/2#max(ys)/10 #y1590 # magnitude
        Gamma2 = 150 # FWHM
        cen2 = 1360
        ming2 =5
        #print(I02)
        
        I03 = y1580*1.1#max(ys)/2#max(ys)/10 #y1590 # magnitude
        Gamma3 = 100 # FWHM
        cen3 = 1590
        q = -6
        
        bnds = ((0.1, I01*5), (cen1-10,cen1+10), (1,10),
                (0.1, I02), (cen2-25,cen2+25), (ming2,300),
                (0.1, I03), (cen3-55,cen3+55), (1,300), #15
                (-100,0), (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'q', 'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]              


    elif material_type == 'Lasered Diamond Doped':     

        y1560 = ys[np.where(xs>1580)[0][-1]]
        print(y1560)
        
        # print('t')
        y2 = np.mean(ys[:20])
        y1 = np.mean(ys[-20:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        # I01 = max(ys) # magnitude
        # Gamma1 = 3 # FWHM
        # cen1 = 1332#1500
        
        # I02 = 1
        # Gamma2 = 2 # FWHM
        # cen2 = 1350
        
        # I03 = 1
        # Gamma3 = 2# FWHM
        # cen3 = 1590
        # q = -6
        
        #if second_pass == True:
            
        I01 = max(ys)#/10 # magnitude
        Gamma1 = 2 # FWHM
        cen1 = 1332#1500
        
        I02 = y1560*1.2##max(ys)/10 #y1590 # magnitude
        Gamma2 = 200 # FWHM
        cen2 = 1300
        
        I03 = y1560*1.2#max(ys)/2#max(ys)/10 #y1590 # magnitude
        Gamma3 = 100 # FWHM
        cen3 = 1555
        q = -6
        
        bnds = ((0.1, I01*2), (cen1-10,cen1+10), (1,10),
                (0.1, I02), (cen2-55,cen2+55), (1,np.inf),
                (0.1, I03), (cen3-55,cen3+55), (1,np.inf), #15
                (-100,0), (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'q', 'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]   
            
            
    elif material_type == 'ie4':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        # y1590 = ys[np.where(xs>1590)[0][-1]] # finds the y value closest to
        #                                      # x = 1590

        I01 = max(ys) # magnitude
        Gamma1 = 4 # FWHM
        cen1 = 1333#1500
        
        I02 = 10#max(ys) # magnitude
        Gamma2 = 10 # FWHM
        cen2 = 1350
        
        I03 = 10 # magnitude
        Gamma3 = 10 # FWHM
        cen3 = 1550
        
        I04 = 10 # magnitude
        Gamma4 = 10 # FWHM
        cen4 = 1590
        q = -6
            
        
        if second_pass == True:
            I01 = max(ys) # magnitude
            Gamma1 = 4 # FWHM
            cen1 = 1333#1500
            
            I02 = max(ys) # magnitude
            Gamma2 = 100 # FWHM
            cen2 = 1350
            
            I03 = max(ys) # magnitude
            Gamma3 = 150 # FWHM
            cen3 = 1550
            
            I04 = max(ys) # magnitude
            Gamma4 = 100 # FWHM
            cen4 = 1590
            q = -6
        
        bnds = ((0.01, I01*2), (cen1-10,cen1+10), (1,10),
                (0.01, I02*5), (cen2-25,cen2+25), (1,250),
                (0.01, I03*5), (cen3-25,cen3+25), (1,200),
                (0.01, I04*5), (cen4-30,cen4+30), (1,250),
                (-100,0), (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, I04, cen4, Gamma4, q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'q', 'm', 'c']  
            
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i] 
            
            
    elif material_type == 'ie5':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        y1560 = ys[np.where(xs>1560)[0][-1]] # finds the y value closest to
                                              # x = 1590
                                              
        print(y1560)
        
        I01 = y1560 # magnitude
        Gamma1 = 20 # FWHM
        cen1 = 1160#1500
        
        I02 = y1560 #y1590 # magnitude
        Gamma2 = 30 # FWHM
        cen2 = 1275
        
        I03 = max(ys) # magnitude
        Gamma3 = 3 # FWHM
        cen3 = 1332#1350
        
        I04 = y1560 # magnitude
        Gamma4 = 30 # FWHM
        cen4 = 1425#1450
        
        I05 = y1560 # magnitude
        Gamma5 = 80# FWHM
        cen5 = 1555
        q = -4
        
        # if second_pass == True:
        #     I01 = max(ys)/10 # magnitude
        #     Gamma1 = 20 # FWHM
        #     cen1 = 1180#1500
            
        #     I02 = max(ys)/10 #y1590 # magnitude
        #     Gamma2 = 30 # FWHM
        #     cen2 = 1310
            
        #     I03 = max(ys)/40 # magnitude
        #     Gamma3 = 1 # FWHM
        #     cen3 = 1850#1500
            
        #     I04 = max(ys)/10 #y1590 # magnitude
        #     Gamma4 = 30 # FWHM
        #     cen4 = 1450
            
        #     I05 = max(ys)/10# y1590 # magnitude
        #     Gamma5 = 30# FWHM
        #     cen5 = 1600
        #     q = -6
        
        bnds = ((0.1, I01*2), (cen1-55,cen1+55), (0,300),
                (0.1, I02*2), (cen2-55,cen2+55), (0,300), 
                (0.1, I03*2), (cen3-5,cen3+5), (0,15),
                (0.1, I04*2), (cen4-55,cen4+55), (0,300),
                (0.1, I05*2), (cen5-55,cen5+55), (0,300), #15
                (-100,0), (m-0.01,m+0.01), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, I04, cen4, Gamma4,
                I05, cen5, Gamma5, q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'q', 'm', 'c']  
            
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i] 
            
            
    elif material_type == 'ie6':     
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)     
        c = float(y2 - m*x2)   
        
        # y1590 = ys[np.where(xs>1590)[0][-1]] # finds the y value closest to
        #                                      # x = 1590
 
        
        I01 = max(ys)/10 # magnitude
        Gamma1 = 20 # FWHM
        cen1 = 1170#1500
        
        I02 = max(ys)/10 #y1590 # magnitude
        Gamma2 = 30 # FWHM
        cen2 = 1260
        
        I03 = max(ys) # magnitude
        Gamma3 = 3 # FWHM
        cen3 = 1340#1500
        
        I04 = max(ys)/10 #y1590 # magnitude
        Gamma4 = 30 # FWHM
        cen4 = 1450

        I05 = max(ys)/10 #y1590 # magnitude
        Gamma5 = 30 # FWHM
        cen5 = 1510
        
        I06 = max(ys)/10# y1590 # magnitude
        Gamma6 = 30# FWHM
        cen6 = 1560
        q = -6
        
        # if second_pass == True:
        #     I01 = max(ys)/10 # magnitude
        #     Gamma1 = 20 # FWHM
        #     cen1 = 1180#1500
            
        #     I02 = max(ys)/10 #y1590 # magnitude
        #     Gamma2 = 30 # FWHM
        #     cen2 = 1310
            
        #     I03 = max(ys)/40 # magnitude
        #     Gamma3 = 1 # FWHM
        #     cen3 = 1850#1500
            
        #     I04 = max(ys)/10 #y1590 # magnitude
        #     Gamma4 = 30 # FWHM
        #     cen4 = 1480
            
        #     I05 = max(ys)/10 #y1590 # magnitude
        #     Gamma5 = 30 # FWHM
        #     cen5 = 1500
            
        #     I06 = max(ys)/10# y1590 # magnitude
        #     Gamma6 = 30# FWHM
        #     cen6 = 1600
        #     q = -6
        
        bnds = ((0.1, I01*5), (cen1-55,cen1+55), (0,np.inf),
                (0.1, I02*5), (cen2-55,cen2+55), (0,np.inf), 
                (0.1, I03*5), (cen3-55,cen3+55), (0,200),
                (0.1, I04*5), (cen4-55,cen4+55), (0,np.inf),
                (0.1, I05*5), (cen5-20,cen5+20), (0,200), #15
                (0.1, I06*5), (cen6-20,cen6+20), (0,200),
                (-100,0), (m-0.1,m+0.1), (0,np.inf))

        # order and shape
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, I04, cen4, Gamma4,
                I05, cen5, Gamma5, I06, cen6, Gamma6, q, m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6', 'q', 'm', 'c']  
            
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i] 

    
    elif material_type == 'IE9':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys)/10 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1200 #1165
        
        I02 = max(ys) # magnitude
        Gamma2 = 100 # FWHM
        cen2 = 1340 #1330
        
        I03 = max(ys)/2 # magnitude
        Gamma3 = 100 # FWHM
        cen3 = 1572 #1505
    
        I04 = max(ys)/2 # magnitude
        Gamma4 = 90 # FWHM
        cen4 = 1590
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (10,200),
                (0, max(ys)*2), (cen2-40,cen2+40), (2,300),
                (0, max(ys)*2), (cen3-50,cen3+50), (2,500),
                (0, max(ys)*2), (cen4-30,cen4+40), (2,300),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]
            
            
            
    if material_type == 'Active Carbon':
    
        #---Lorentzian D2 peak---
        lowerx = np.where(xs>1100)
        upperx = np.where(xs<1180)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D2 peak
        Gamma1 = 40 # half width at half maximum (HWHM)
        cen1 = 1167 # center of the D2 peak   
        
        #---Lorentzian D peak---
        lowerx = np.where(xs>1250)
        upperx = np.where(xs<1400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I02 = max(ys[interx]) # magnitude of the D peak
        Gamma2 = 55 # half width at half maximum (HWHM)
        cen2 = 1315 # center of the D peak
        
        #---BWF G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1630)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I03 = max(ys[interx])  # magnitude of the G peak
        Gamma3 = 70 # half width at half maximum (HWHM)
        cen3 = 1590 # center of the G peak
        q = -4 # where 1/q is the Fano parameter
        

        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        #---Optimization Bounds---
        bnds = ((0,I01), (cen1-30,cen1+30), (0,np.inf), 
                (0, I02+200), (cen2-30,cen2+50), (0,np.inf),
                (0,I03+200), (cen3-50, cen3+30), (0,np.inf),
                (-500,-0.1), (m-0.1,m+0.1), (0,np.inf)) # happy w/ m+-0.1 as this should capture the extremes in the noise
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)   
        
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, q, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'q', 'm', 'c']
        
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    if material_type == 'AC kc8':
    
        #---Lorentzian D2 peak---
        lowerx = np.where(xs>1100)
        upperx = np.where(xs<1180)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D2 peak
        Gamma1 = 40 # half width at half maximum (HWHM)
        cen1 = 1600#1167 # center of the D2 peak   
        
        #---Lorentzian D peak---
        lowerx = np.where(xs>1250)
        upperx = np.where(xs<1400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I02 = max(ys[interx]) # magnitude of the D peak
        Gamma2 = 55 # half width at half maximum (HWHM)
        cen2 = 1315 # center of the D peak
        
        #---BWF G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1630)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I03 = max(ys[interx])  # magnitude of the G peak
        Gamma3 = 60 # half width at half maximum (HWHM)
        cen3 = 1540 # center of the G peak
        q = -40
        
        #---Lorentzian G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1630)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I04 = max(ys[interx])  # magnitude of the G peak
        Gamma4 = 60 # half width at half maximum (HWHM)
        cen4 = 1580 # center of the G peak
        
        # if second_pass == True:
        #     #---BWF G peak---
        #     lowerx = np.where(xs>1570)
        #     upperx = np.where(xs<1630)
        #     interx = np.intersect1d(lowerx,upperx) # indices for desired range
        #     I03 = max(ys[interx])  # magnitude of the G peak
        #     Gamma3 = 60 # half width at half maximum (HWHM)
        #     cen3 = 1580 # center of the G peak
        #     q = -40
            
        #     #---Lorentzian G peak---
        #     lowerx = np.where(xs>1570)
        #     upperx = np.where(xs<1630)
        #     interx = np.intersect1d(lowerx,upperx) # indices for desired range
        #     I04 = 5#max(ys[interx])  # magnitude of the G peak
        #     Gamma4 = 3 # half width at half maximum (HWHM)
        #     cen4 = 1950 # center of the G peak            
                
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        #---Optimization Bounds---
        bnds = ((0,I01), (cen1-30,cen1+30), (0,np.inf), 
                (0, I02+200), (cen2-30,cen2+30), (0,np.inf),
                (0,I03+200), (cen3-20, cen3+20), (0,np.inf),
                (0, I04+200), (cen4-30, cen4+30), (0,np.inf),
                (-500,-0.1), (m-0.01,m+0.01), (0,np.inf)) # happy w/ m+-0.1 as this should capture the extremes in the noise
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)   
        
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, q, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'q', 'm', 'c']
        
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    if material_type == 'AC MAD mix':
    
        #---Lorentzian D2 peak---
        lowerx = np.where(xs>1100)
        upperx = np.where(xs<1180)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D2 peak
        Gamma1 = 40 # half width at half maximum (HWHM)
        cen1 = 1167 # center of the D2 peak   
        
        #---Lorentzian D peak---
        lowerx = np.where(xs>1250)
        upperx = np.where(xs<1400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I02 = max(ys[interx]) # magnitude of the D peak
        Gamma2 = 55 # half width at half maximum (HWHM)
        cen2 = 1315 # center of the D peak
        
        #---BWF G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1630)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I03 = max(ys[interx])  # magnitude of the G peak
        Gamma3 = 60 # half width at half maximum (HWHM)
        cen3 = 1540 # center of the G peak
        q = -40
        
        #---Lorentzian G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1630)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I04 = max(ys[interx])  # magnitude of the G peak
        Gamma4 = 60 # half width at half maximum (HWHM)
        cen4 = 1580 # center of the G peak
                
        #---Lorentzian G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1630)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I05 = max(ys[interx])  # magnitude of the G peak
        Gamma5 = 60 # half width at half maximum (HWHM)
        cen5 = 1600 # center of the G peak               
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        #---Optimization Bounds---
        bnds = ((0,I01), (cen1-30,cen1+30), (0,np.inf), 
                (0, I02+200), (cen2-30,cen2+30), (0,np.inf),
                (0,I03+200), (cen3-20, cen3+20), (0,np.inf),
                (0, I04+200), (cen4-30, cen4+30), (0,np.inf),
                (0, I05+200), (cen5-30, cen5+30), (0,np.inf),
                (-500,-0.1), (m-0.01,m+0.01), (0,np.inf)) # happy w/ m+-0.1 as this should capture the extremes in the noise
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)   
        
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, q, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4',
               'I05', 'cen5', 'Gamma5', 'q', 'm', 'c']
        
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]


    if material_type == 'Carbon Black 3p':
      
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)        

        #---Lorentzian D peak---
        lowerx = np.where(xs>1250)
        upperx = np.where(xs<1400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D peak
        Gamma1 = 55 # half width at half maximum (HWHM)
        cen1 = 1315 # center of the D peak
        
        #---Lorentzian D2 peak---
        lowerx = np.where(xs>1100)
        upperx = np.where(xs<1180)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I02 = max(ys[interx]) # magnitude of the D2 peak
        Gamma2 = 40 # half width at half maximum (HWHM)
        cen2 = 1500 # center of the D2 peak   
        
        #---BWF G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1630)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I03 = max(ys[interx])  # magnitude of the G peak
        Gamma3 = 25 # half width at half maximum (HWHM)
        cen3 = 1590 # center of the G peak
        q = -4 # where 1/q is the Fano parameter
        
        
        #---Optimization Bounds---
        bnds = ((0, I01+200), (cen1-30,cen1+60), (0,np.inf),
                (0, I02), (cen2-30,cen2+30), (0,250), 
                (0 ,I03+200), (cen3-30, cen3+30), (0,np.inf),
                (-500,-0.1), (m-0.01,m+0.01), (0,np.inf)) # happy w/ m+-0.1 as this should capture the extremes in the noise
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)   
        
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, q, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'q', 'm', 'c']
        
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]        
            
    if material_type == 'Carbon Black 4p':
      
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)        

        #---Lorentzian D2 peak---
        I01 = max(ys)/5 # magnitude of the D2 peak
        Gamma1 = 60 # half width at half maximum (HWHM)
        cen1 = 1170#67 # center of the D2 peak   

        #---Lorentzian D peak---
        I02 = max(ys) # magnitude of the D peak
        Gamma2 = 55 # half width at half maximum (HWHM)
        cen2 = 1315 # center of the D peak
        
        #---BWF G peak---
        I03 = max(ys)  # magnitude of the G peak
        Gamma3 = 55 # half width at half maximum (HWHM)
        cen3 = 1600 # center of the G peak
        q = -4 # where 1/q is the Fano parameter
        
        #---Lorentzian A peak---
        I04 = max(ys)/4 # magnitude of the D peak
        Gamma4 = 55 # half width at half maximum (HWHM)
        cen4 = 1490 # center of the D peak
        
        
        #---Optimization Bounds---
        bnds = ((0, I02*2), (cen1-40,cen1+40), (0,np.inf),
                (0, I02*2), (cen2-30,cen2+60), (0,np.inf), 
                (0 ,I02*2), (cen3-20, cen3+20), (10,np.inf),
                (0 ,I02*2), (cen4-25, cen4+25), (0,250),
                (-500,-0.1), (m-0.01,m+0.01), (0,np.inf)) # happy w/ m+-0.1 as this should capture the extremes in the noise
        bnds = np.array(bnds)
        bnds = np.transpose(bnds)   
        
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, q, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4',
               'q', 'm', 'c']
        
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]     
            
        
    elif material_type == 'MWCNT':
        
        #---Lorentzian D peak---
        lowerx = np.where(xs>1320)
        upperx = np.where(xs<1400)
        interx = np.intersect1d(lowerx,upperx)
        I0D = max(ys[interx]) # magnitude of the D peak
        gammaD = 40 # half width at half maximum (HWHM)
        cenD = 1360
        
        #---Lorentzian G peak---
        lowerx = np.where(xs>1570)
        upperx = np.where(xs<1610)
        interx = np.intersect1d(lowerx,upperx)
        I0G = max(ys[interx])  # magnitude of the G peak
        gammaG = 40 # half width at half maximum (HWHM)
        cenG = 1590
                    
        #---Background---
        lowerx = np.where(xs<1250)
        upperx = np.where(xs>1750)
        y2 = np.mean(ys[upperx])
        y1 = np.mean(ys[lowerx])
        x2 = max(xs)
        x1 = min(xs)
        m = float((y2-y1)/(x2-x1))
        c = float(y2 - m*x2)
        
        #---Optimization Bounds---
        bnds = ((0, I0D+400), (cenD-20,cenD+20), (0,np.inf),
                (0,I0G+400), (cenG-20, cenG+20), (0,np.inf), 
                (-np.inf,np.inf), (-np.inf,np.inf))
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        guesses = [I0D, cenD, gammaD,
                I0G, cenG, gammaG,
                m, c]
        p_s = ['I0D', 'cenD', 'gammaD', 
               'I0G', 'cenG', 'gammaG', 
               'm', 'c']

        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'CD 785 4L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1330
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1480 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1590  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1170 #1560-1630
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (2,300),
                (0, max(ys)*2), (cen2-60,cen2+60), (2,300),
                (0, max(ys)*2), (cen3-50,cen3+80), (2,500),
                (0, max(ys)*2), (cen4-30,cen4+40), (10,300),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]


    elif material_type == 'CD 785 5L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1170
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1300 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1585  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1610 #1560-1630
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 2650 #1560-1630
        
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (50,300),
                (0, max(ys)*2), (cen2-60,cen2+60), (2,300),
                (0, max(ys)*2), (cen3-10,cen3+10), (2,500),
                (0, max(ys)*2), (cen4-10,cen4+10), (2,300),
                (0, max(ys)*2), (cen5-50,cen5+80), (2,500),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'm', 'c',]
        
        
    elif material_type == 'CD 785 7L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1300
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1450 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1585  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1610 #1560-1630
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 2650 #1560-1630
        
        I06 = max(ys)/2 # magnitude
        Gamma6 = 80 # FWHM
        cen6 = 2850 #1560-1630
        
        I07 = max(ys)/2 # magnitude
        Gamma7 = 80 # FWHM
        cen7 = 2950 #1560-1630
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (50,300),
                (0, max(ys)*2), (cen2-60,cen2+60), (2,300),
                (0, max(ys)*2), (cen3-50,cen3+80), (2,500),
                (0, max(ys)*2), (cen4-30,cen4+40), (2,300),
                (0, max(ys)*2), (cen5-50,cen5+80), (2,500),
                (0, max(ys)*2), (cen6-30,cen6+40), (2,300),
                (0, max(ys)*2), (cen7-50,cen7+80), (2,500),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7, m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
               'I07', 'cen7', 'Gamma7', 'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'CD 514 7L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1355
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1450

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1585  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1630 #1560-1630
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 2720 #1560-1630
        
        I06 = max(ys)/2 # magnitude
        Gamma6 = 80 # FWHM
        cen6 = 2850 #1560-1630
        
        I07 = max(ys)/2 # magnitude
        Gamma7 = 80 # FWHM
        cen7 = 2950 #1560-1630
        
        bnds = ((0, max(ys)*2), (cen1-30,cen1+30), (2,300),
                (0, max(ys)*2), (cen2-50,cen2+40), (2,300),
                (0, max(ys)*2), (cen3-10,cen3+10), (2,500),
                (0, max(ys)*2), (cen4-30,cen4+30), (2,300),
                (0, max(ys)*2), (cen5-50,cen5+80), (2,500),
                (0, max(ys)*2), (cen6-30,cen6+40), (2,300),
                (0, max(ys)*2), (cen7-20,cen7+20), (2,500),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7, m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
               'I07', 'cen7', 'Gamma7', 'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]
            
            
    elif material_type == 'CD 514 8L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1330
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1355

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1450  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1585 #1560-1630
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 1615 #1560-1630
        
        I06 = max(ys)/2 # magnitude
        Gamma6 = 80 # FWHM
        cen6 = 2720 #1560-1630
        
        I07 = max(ys)/2 # magnitude
        Gamma7 = 80 # FWHM
        cen7 = 2850 #1560-1630
        
        I08 = max(ys)/2 # magnitude
        Gamma8 = 80 # FWHM
        cen8 = 2950 #1560-1630
        
        bnds = ((0, max(ys)*2), (cen1-10,cen1+10), (10,300),
                (0, max(ys)*2), (cen2-10,cen2+10), (2,300),
                (0, max(ys)*2), (cen3-10,cen3+10), (2,500),
                (0, max(ys)*2), (cen4-10,cen4+10), (2,300),
                (0, max(ys)*2), (cen5-50,cen5+80), (2,500),
                (0, max(ys)*2), (cen6-30,cen6+40), (2,300),
                (0, max(ys)*2), (cen7-20,cen7+20), (2,500),
                (0, max(ys)*2), (cen8-20,cen8+20), (2,300),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7, I08, cen8, Gamma8,
                   m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
               'I07', 'cen7', 'Gamma7', 'I08', 'cen8', 'Gamma8',
               'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]
            
    elif material_type == 'CD 785 8L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1300
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1445 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1460  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1585 #1560-1630
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 1610 #1560-1630
        
        I06 = max(ys)/2 # magnitude
        Gamma6 = 80 # FWHM
        cen6 = 2650 #1560-1630
        
        I07 = max(ys)/2 # magnitude
        Gamma7 = 80 # FWHM
        cen7 = 2850 #1560-1630
        
        I08 = max(ys)/2 # magnitude
        Gamma8 = 80 # FWHM
        cen8 = 2950 #1560-1630
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (50,300),
                (0, max(ys)*2), (cen2-60,cen2+60), (2,300),
                (0, max(ys)*2), (cen3-50,cen3+80), (2,500),
                (0, max(ys)*2), (cen4-30,cen4+40), (2,300),
                (0, max(ys)*2), (cen5-50,cen5+80), (2,500),
                (0, max(ys)*2), (cen6-30,cen6+40), (2,300),
                (0, max(ys)*2), (cen7-50,cen7+80), (2,500),
                (0, max(ys)*2), (cen8-30,cen8+40), (2,300),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7, I08, cen8, Gamma8,
                   m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
               'I07', 'cen7', 'Gamma7', 'I08', 'cen8', 'Gamma8',
               'm', 'c',]
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'CD 514 4L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys)/10 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1280 #1165
        
        I02 = max(ys) # magnitude
        Gamma2 = 100 # FWHM
        cen2 = 1370 #1330
        
        I03 = max(ys)/2 # magnitude
        Gamma3 = 100 # FWHM
        cen3 = 1520 #1505
    
        I04 = max(ys)/2 # magnitude
        Gamma4 = 90 # FWHM
        cen4 = 1590
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (10,200),
                (0, max(ys)*2), (cen2-40,cen2+40), (2,300),
                (0, max(ys)*2), (cen3-50,cen3+50), (2,500),
                (0, max(ys)*2), (cen4-30,cen4+40), (2,300),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'CD 514 5L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1240
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1370 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1520  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1590 #1560-1630
        
        I05 = max(ys)/20 # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 1620 #1560-1630
        
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+20), (2,200),
                (0, max(ys)*2), (cen2-30,cen2+30), (2,300),
                (0, max(ys)*2), (cen3-40,cen3+40), (2,500),
                (0, max(ys)*2), (cen4-15,cen4+15), (2,300),
                (0, max(ys)*2), (cen5-15,cen5+15), (2,300),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4',
                'I05', 'cen5', 'Gamma5', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]
            
    elif material_type == 'CD 514 6L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1360
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1585 

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1610  #1470-1570 # 1430-1520 #1430-1560, so 1495 pm 65
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 2250 #1560-1630
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 2690 #1560-1630
        
        I06 = max(ys)/2 # magnitude
        Gamma6 = 80 # FWHM
        cen6 = 2950 #1560-1630
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (2,300),
                (0, max(ys)*2), (cen2-10,cen2+10), (2,300),
                (0, max(ys)*2), (cen3-10,cen3+20), (2,500),
                (0, max(ys)*2), (cen4-40,cen4+40), (2,300),
                (0, max(ys)*2), (cen5-40,cen5+40), (2,300),
                (0, max(ys)*2), (cen6-40,cen6+40), (2,300),
            (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4',
                'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]

    elif material_type == 'CD Grain 7 3L':
        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1330

        I02 = max(ys)/10 # magnitude
        Gamma2 = 75 # FWHM
        cen2 = 1480   
        # # I02 = 1.5 # magnitude low (as if peak not present)
        # # Gamma2 = 2 # FWHM
        # # cen2 = 1000
        # if second_pass == True:
        #     I02 = max(ys)/10 # magnitude
        #     Gamma2 = 75 # FWHM
        #     cen2 = 1500            
        
        I03 = max(ys)/2 # magnitude
        Gamma3 = 80 # FWHM
        cen3 = 1590
        
        bnds = ((0, max(ys)*2), (cen1-30,cen1+40), (2,200),
                (0, max(ys)*2), (cen2-70,cen2+40), (2,500),
                (0, max(ys)*2), (cen3-70,cen3+20), (2,300),
            (m-0.001,m+0.001), (0,np.inf)) # g peak lims +-20 
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                I03, cen3, Gamma3, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]
            
    elif material_type == 'CD Grain 7 4L':
          
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1165
        
        I02 = max(ys) # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1330

        I03 = max(ys)/10 # magnitude
        Gamma3 = 75 # FWHM
        cen3 = 1480   
        # if second_pass == True:
        #     I03 = max(ys)/10 # magnitude
        #     Gamma3 = 75 # FWHM
        #     cen3 = 1490            
        
        I04 = max(ys)/2 # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1590
        
        bnds = ((0, max(ys)*2), (cen1-60,cen1+60), (10,200),
                (0, max(ys)*2), (cen2-30,cen2+40), (2,200),
                (0, max(ys)*2), (cen3-70,cen3+40), (2,500),
                (0, max(ys)*2), (cen4-30,cen4+20), (2,300),
            (m-0.001,m+0.001), (0,np.inf)) # g peak lims +-20 
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]


    elif material_type == 'CD Grain 6 2L':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        I01 = max(ys) # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1380
        
        I02 = max(ys)/2 # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1567
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+50), (2,600),
                (0, max(ys)*2), (cen2-50,cen2+60), (2,600),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            
    elif material_type == 'CD Grain 6 3L':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        if sample == 'CD G6 514':
            I01 = max(ys)/3 # magnitude
            Gamma1 = 80 # FWHM
            cen1 = 1376#1165
            
            I02 = max(ys) # magnitude
            Gamma2 = 80 # FWHM
            cen2 = 1570#1330
            
            I03 = max(ys)/2 # magnitude
            Gamma3 = 80 # FWHM
            cen3 = 1590#1550
            
        elif sample == 'CD G6 785':
            I01 = max(ys)/3 # magnitude
            Gamma1 = 80 # FWHM
            cen1 = 1330#1165
            
            I02 = max(ys) # magnitude
            Gamma2 = 80 # FWHM
            cen2 = 1505#1330
            
            I03 = max(ys)/2 # magnitude
            Gamma3 = 80 # FWHM
            cen3 = 1590#1550
        
        bnds = ((0, max(ys)), (cen1-50,cen1+50), (2,600),
                (0, max(ys)*2), (cen2-100,cen2+100), (2,600),
                (0, max(ys)*2), (cen3-100,cen3+100), (2,600),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            
    elif material_type == 'CD Grain 6 4L':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        if sample == 'CD G6 785':
            I01 = max(ys)/20 # magnitude
            Gamma1 = 40 # FWHM
            cen1 = 1165 #1165
            
            I02 = max(ys) # magnitude
            Gamma2 = 100 # FWHM
            cen2 = 1330 #1330
            
            I03 = max(ys)/2 # magnitude
            Gamma3 = 100 # FWHM
            cen3 = 1505 #1505
        
            I04 = max(ys)/2 # magnitude
            Gamma4 = 90 # FWHM
            cen4 = 1590
 

            
        elif sample == 'CD G6 514':
            I01 = max(ys)/90 # magnitude
            Gamma1 = 4 # FWHM
            cen1 = 1240 #1165
            
            I02 = max(ys) # magnitude
            Gamma2 = 100 # FWHM
            cen2 = 1370 #1330
            
            I03 = max(ys)/2 # magnitude
            Gamma3 = 100 # FWHM
            cen3 = 1520 #1505
        
            I04 = max(ys)/2 # magnitude
            Gamma4 = 90 # FWHM
            cen4 = 1590
        
        bnds = ((0, max(ys)*2), (cen1-50,cen1+10), (2,200),
                (0, max(ys)*2), (cen2-30,cen2+30), (2,300),
                (0, max(ys)*2), (cen3-50,cen3+50), (2,400),
                (0, max(ys)*2), (cen4-50,cen4+50), (2,120),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            
    elif material_type == 'CD Grain 6 5L':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys)/3 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1370 #1165
        
        I02 = max(ys) # magnitude
        Gamma2 = 100 # FWHM
        cen2 = 1510 #1330
        
        I03 = max(ys)/2 # magnitude
        Gamma3 = 100 # FWHM
        cen3 = 1565 #1505

        I04 = max(ys)/2 # magnitude
        Gamma4 = 90 #50 # FWHM
        cen4 = 1590
        
        I05 = max(ys)/2 # magnitude
        Gamma5 = 90 #50 # FWHM
        cen5 = 1610
        

        
        bnds = ((0, max(ys)), (cen1-100,cen1+100), (2,600),
                (0, max(ys)*2), (cen2-50,cen2+50), (2,600),
                (0, max(ys)*2), (cen3-100,cen3+100), (2,600),
                (0, max(ys)*2), (cen4-100,cen4+100), (2,200),
                (0, max(ys)*2), (cen5-100,cen5+100), (2,200),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                'I05', 'cen5', 'Gamma5', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    

    elif material_type == 'CD Grain 6 2V':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        amp1 = max(ys)/6
        sigma1 = 10
        I01 = max(ys)/2 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1376
    
        amp2 = max(ys)/6
        sigma2 = 10
        I02 = max(ys)/2 # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1590

        bnds = ((0, np.inf), (0, np.inf), (0, np.inf), (cen1-50,cen1+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen2-50,cen2+50), (0,np.inf),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 

        test = dict()
        
        guesses = [amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2, m, c]

        p_s = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
               'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2', 'm', 'c']
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            
    elif material_type == 'CD Grain 6 3V':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        amp1 = max(ys)/6
        sigma1 = 10
        I01 = max(ys)/2 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1376
    
        amp2 = max(ys)/6
        sigma2 = 10
        I02 = max(ys)/2 # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1165
        
        amp3 = max(ys)/6
        sigma3 = 10
        I03 = max(ys) # magnitude
        Gamma3 = 80 # FWHM
        cen3 = 1590

        bnds = ((0, np.inf), (0, np.inf), (0, np.inf), (cen1-50,cen1+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen2-50,cen2+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen3-50,cen3+50), (0,np.inf),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 

        test = dict()
        
        guesses = [amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2,
                   amp3, sigma3, I03, cen3, Gamma3, m, c]

        p_s = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
               'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2',
               'amp3', 'sigma3', 'I03', 'cen3', 'Gamma3', 'm', 'c']
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            
            
    elif material_type == 'CD Grain 6 4V':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        amp1 = max(ys)/6
        sigma1 = 10
        I01 = max(ys)/2 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1165
    
        amp2 = max(ys)/6
        sigma2 = 10
        I02 = max(ys)/2 # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1375
        
        amp3 = max(ys)/6
        sigma3 = 10
        I03 = max(ys) # magnitude
        Gamma3 = 80 # FWHM
        cen3 = 1550
        
        amp4 = max(ys)/6
        sigma4 = 10
        I04 = max(ys) # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1590

        bnds = ((0, np.inf), (0, np.inf), (0, np.inf), (cen1-50,cen1+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen2-50,cen2+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen3-50,cen3+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen4-50,cen4+50), (0,np.inf),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 

        test = dict()
        
        guesses = [amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2,
                   amp3, sigma3, I03, cen3, Gamma3,
                   amp4, sigma4, I04, cen4, Gamma4, m, c]

        p_s = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
               'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2',
               'amp3', 'sigma3', 'I03', 'cen3', 'Gamma3',
               'amp4', 'sigma4', 'I04', 'cen4', 'Gamma4', 'm', 'c']
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    

    elif material_type == 'CD Grain 6 5V':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        amp1 = max(ys)
        sigma1 = 1000
        I01 = max(ys)/2 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1165
    
        amp2 = max(ys)
        sigma2 = 1000
        I02 = max(ys)/2 # magnitude
        Gamma2 = 80 # FWHM
        cen2 = 1375
        
        amp3 = max(ys)
        sigma3 = 1000
        I03 = max(ys) # magnitude
        Gamma3 = 80 # FWHM
        cen3 = 1550
        
        amp4 = max(ys)
        sigma4 = 1000
        I04 = max(ys) # magnitude
        Gamma4 = 80 # FWHM
        cen4 = 1590
        
        amp5 = max(ys)
        sigma5 = 1000
        I05 = max(ys) # magnitude
        Gamma5 = 80 # FWHM
        cen5 = 1610

        bnds = ((0, np.inf), (0, np.inf), (0, np.inf), (cen1-50,cen1+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen2-50,cen2+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen3-50,cen3+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen4-50,cen4+50), (0,np.inf),
                (0, np.inf), (0, np.inf), (0, np.inf), (cen5-50,cen5+50), (0,np.inf),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 

        test = dict()
        
        guesses = [amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2,
                   amp3, sigma3, I03, cen3, Gamma3,
                   amp4, sigma4, I04, cen4, Gamma4, 
                   amp5, sigma5, I05, cen5, Gamma5, m, c]

        p_s = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
               'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2',
               'amp3', 'sigma3', 'I03', 'cen3', 'Gamma3',
               'amp4', 'sigma4', 'I04', 'cen4', 'Gamma4',
               'amp5', 'sigma5', 'I05', 'cen5', 'Gamma5', 'm', 'c']
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    

    elif material_type == 'CD Grain 8 2L':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys)/3 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1320
        
        I02 = max(ys) # magnitude
        Gamma2 = 3 # FWHM
        cen2 = 1600

        bnds = ((0, max(ys)*2), (cen1-60,cen1+60), (2,300),
                (0, max(ys)*2), (cen2-50,cen2+50), (2,200),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            

    elif material_type == 'CD Grain 8 3L':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys)/3 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1320
        
        I02 = max(ys) # magnitude
        Gamma2 = 3 # FWHM
        cen2 = 1585
        
        I03 = max(ys)/2 # magnitude
        Gamma3 = 100 # FWHM
        cen3 = 2700#1610

        bnds = ((0, max(ys)*2), (cen1-60,cen1+60), (2,300),
                (0, max(ys)*2), (cen2-50,cen2+50), (2,200),
                (0, max(ys)*2), (cen3-50,cen3+50), (2,400),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            
            
    elif material_type == 'CD Grain 8 4L':
        
        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)

        I01 = max(ys)/3 # magnitude
        Gamma1 = 80 # FWHM
        cen1 = 1320
        
        I02 = max(ys) # magnitude
        Gamma2 = 3 # FWHM
        cen2 = 1333
        
        I03 = max(ys)/2 # magnitude
        Gamma3 = 100 # FWHM
        cen3 = 1585

        I04 = max(ys)/4 # magnitude
        Gamma4 = 30 # FWHM
        cen4 = 1610

        bnds = ((0, max(ys)*2), (cen1-60,cen1+60), (2,300),
                (0, max(ys)*2), (cen2-50,cen2+50), (2,10),
                (0, max(ys)*2), (cen3-50,cen3+50), (2,400),
                (0, max(ys)*2), (cen4-50,cen4+50), (2,120),
                (m-0.001,m+0.001), (0,np.inf))
        
        bnds = np.array(bnds)
        bnds = np.transpose(bnds) 
        
        test = dict()
        
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   m, c]    
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4',
               'm', 'c']
    
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]    
            

    elif material_type == 'MnOx':
        
        #---Lorentzian 575 peak---
        lowerx = np.where(xs>500)
        upperx = np.where(xs<600)
        interx = np.intersect1d(lowerx,upperx)
        I03 = max(ys[interx]) # magnitude of the D peak
        Gamma3 = 20 # half width at half maximum (HWHM)
        cen3 = 565
        
        #---Lorentzian 635 peak---
        lowerx = np.where(xs>600)
        upperx = np.where(xs<700)
        interx = np.intersect1d(lowerx,upperx)
        I04 = max(ys[interx])  # magnitude of the G peak
        Gamma4 = 20 # half width at half maximum (HWHM)
        cen4 = 631
        
        # #---Lorentzian 541 shoulder---
        # lowerx = np.where(xs>570)
        # upperx = np.where(xs<590)
        # interx = np.intersect1d(lowerx,upperx) # indices for desired range
        # I0D2 = max(ys[interx]) # magnitude of the D2 peak
        # gammaD2 = 10 # half width at half maximum (HWHM)
        # cenD2 = 541 # center of the D2 peak
        
        #---Lorentzian 375 peak---
        lowerx = np.where(xs>300)
        upperx = np.where(xs<400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D2 peak
        Gamma1 = 20 # half width at half maximum (HWHM)
        cen1 = 375 # center of the D2 peak
        
        #---Lorentzian 443 peak---
        I02 = I01 # magnitude of the D2 peak
        Gamma2 = 25 # half width at half maximum (HWHM)
        cen2 = 453 # center of the D2 peak
                    
        #---Background---        
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        #---Optimization Bounds---
        bnds = ((0, I01+2000), (cen1-20,cen1+20), (0,np.inf),
                (0, I02+2000), (cen2-20,cen2+20), (0,np.inf),
                (0, I03+2000), (cen3-10,cen3+10), (0,np.inf),
                (0, I04+2000), (cen4-20, cen4+20), (0,np.inf),
                (m-0.1,m+0.1), (0,np.inf))

        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4, 
                   m, c]
        
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'm', 'c']  
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]
            
    elif material_type == 'MnOx2':

        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        #---Lorentzian 375 peak---
        lowerx = np.where(xs>300)
        upperx = np.where(xs<400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D2 peak
        Gamma1 = 20 # half width at half maximum (HWHM)
        cen1 = 360 # center of the D2 peak
        
        #---Lorentzian 443 peak---
        I02 = I01 # magnitude of the D2 peak
        Gamma2 = 25 # half width at half maximum (HWHM)
        cen2 = 455 # center of the D2 peak
        
        #---Lorentzian 575 peak---
        lowerx = np.where(xs>500)
        upperx = np.where(xs<600)
        interx = np.intersect1d(lowerx,upperx)
        I03 = max(ys[interx]) # magnitude of the D peak
        Gamma3 = 20 # half width at half maximum (HWHM)
        cen3 = 565
        
        #---Lorentzian 635 peak---
        lowerx = np.where(xs>600)
        upperx = np.where(xs<700)
        interx = np.intersect1d(lowerx,upperx)
        I04 = max(ys[interx])  # magnitude of the G peak
        Gamma4 = 20 # half width at half maximum (HWHM)
        cen4 = 625
        
        #---Lorentzian 635 peak---
        lowerx = np.where(xs>600)
        upperx = np.where(xs<700)
        interx = np.intersect1d(lowerx,upperx)
        I05 = max(ys[interx])  # magnitude of the G peak
        I0max = I05+2000
        Gamma5 = 15 # half width at half maximum (HWHM)
        Gmin = 10
        Gmax = 100
        cen5 = 655
        
        if second_pass == True:
            Gamma5 = 3
            Gmin = 2
            Gmax = 10
            I05 = 1
            I0max = 10
        
        #---Optimization Bounds---
        bnds = ((0, I01+2000), (cen1-30,cen1+30), (0,300),
                (0, I02+2000), (cen2-50,cen2+30), (0,300),
                (0, I03+2000), (cen3-15,cen3+15), (0,300),
                (0, I04+2000), (cen4-25, cen4+15), (0,300),
                (0, I0max), (cen5-15, cen5+15), (Gmin,Gmax),
                (m-0.01,m+0.01), (0,np.inf))

        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'm', 'c',] 
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]


    elif material_type == 'MnOx3':

        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        #---Lorentzian 375 peak---
        lowerx = np.where(xs>300)
        upperx = np.where(xs<400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D2 peak
        Gamma1 = 40 # half width at half maximum (HWHM)
        cen1 = 300 # center of the D2 peak
        
        #---Lorentzian 443 peak---
        I02 = I01 # magnitude of the D2 peak
        Gamma2 = 25 # half width at half maximum (HWHM)
        cen2 = 360 # center of the D2 peak
        
        #---Lorentzian 443 peak---
        I03 = I01 # magnitude of the D2 peak
        Gamma3 = 25 # half width at half maximum (HWHM)
        cen3 = 410 # center of the D2 peak
        
        #---Lorentzian 443 peak---
        I04 = I01 # magnitude of the D2 peak
        Gamma4 = 25 # half width at half maximum (HWHM)
        cen4 = 460 # center of the D2 peak
        
        #---Lorentzian 575 peak---
        lowerx = np.where(xs>500)
        upperx = np.where(xs<600)
        interx = np.intersect1d(lowerx,upperx)
        I05 = max(ys[interx]) # magnitude of the D peak
        Gamma5 = 20 # half width at half maximum (HWHM)
        cen5 = 565
        
        #---Lorentzian 635 peak---
        lowerx = np.where(xs>600)
        upperx = np.where(xs<700)
        interx = np.intersect1d(lowerx,upperx)
        I06 = max(ys[interx])  # magnitude of the G peak
        Gamma6 = 20 # half width at half maximum (HWHM)
        cen6 = 625
        
        #---Lorentzian 635 peak---
        lowerx = np.where(xs>600)
        upperx = np.where(xs<700)
        interx = np.intersect1d(lowerx,upperx)
        I07 = max(ys[interx])  # magnitude of the G peak
        I0max = I05+2000
        Gamma7 = 15 # half width at half maximum (HWHM)
        Gmin = 10
        Gmax = 100
        cen7 = 655
        
        if second_pass == True:
            Gamma7 = 3
            Gmin = 2
            Gmax = 10
            I07 = 1
            I0max = 10
        
        #---Optimization Bounds---
        bnds = ((0, I01+2000), (cen1-30,cen1+30), (10,300),
                (0, I02+2000), (cen2-20,cen2+20), (10,300),
                (0, I03+2000), (cen3-20,cen3+20), (10,300),
                (0, I04+2000), (cen4-30,cen4+30), (10,300),
                (0, I05+2000), (cen5-15,cen5+15), (10,300),
                (0, I06+2000), (cen6-25, cen6+15), (10,300),
                (0, I0max), (cen7-15, cen7+15), (Gmin,Gmax),
                (m-0.01,m+0.01), (0,np.inf))

        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7,
                   m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
               'I07', 'cen7', 'Gamma7',
               'm', 'c',]
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]            
 

    elif material_type == 'MnOx4':

        #---Background---
        y2 = np.mean(ys[:10])
        y1 = np.mean(ys[-10:])
        x2 = xs[0]
        x1 = xs[-1]
        m = (y2-y1)/(x2-x1)
        c = float(y2 - m*x2)
        
        #---Lorentzian 375 peak---
        lowerx = np.where(xs>300)
        upperx = np.where(xs<400)
        interx = np.intersect1d(lowerx,upperx) # indices for desired range
        I01 = max(ys[interx]) # magnitude of the D2 peak
        Gamma1 = 40 # half width at half maximum (HWHM)
        cen1 = 340 # center of the D2 peak
        
        #---Lorentzian 443 peak---
        I02 = I01 # magnitude of the D2 peak
        Gamma2 = 25 # half width at half maximum (HWHM)
        cen2 = 390 # center of the D2 peak
        
        #---Lorentzian 443 peak---
        I03 = I01 # magnitude of the D2 peak
        Gamma3 = 25 # half width at half maximum (HWHM)
        cen3 = 460 # center of the D2 peak
        
        #---Lorentzian 575 peak---
        lowerx = np.where(xs>500)
        upperx = np.where(xs<600)
        interx = np.intersect1d(lowerx,upperx)
        I04 = max(ys[interx]) # magnitude of the D peak
        Gamma4 = 20 # half width at half maximum (HWHM)
        cen4 = 565
        
        #---Lorentzian 635 peak---
        lowerx = np.where(xs>600)
        upperx = np.where(xs<700)
        interx = np.intersect1d(lowerx,upperx)
        I05 = max(ys[interx])  # magnitude of the G peak
        Gamma5 = 20 # half width at half maximum (HWHM)
        cen5 = 625
        
        #---Lorentzian 635 peak---
        lowerx = np.where(xs>600)
        upperx = np.where(xs<700)
        interx = np.intersect1d(lowerx,upperx)
        I06 = 20#max(ys[interx])  # magnitude of the G peak
        I0max = I05+2000
        Gamma6 = 15 # half width at half maximum (HWHM)
        Gmin = 10
        Gmax = 100
        cen6 = 655#655
        
        # if second_pass == True:
        #     Gamma6 = 3
        #     Gmin = 2
        #     Gmax = 10
        #     I06 = 1
        #     I0max = 10
        
        #---Optimization Bounds---
        bnds = ((0, I01+2000), (cen1-30,cen1+30), (10,120),
                (0, I02+2000), (cen2-20,cen2+20), (10,120),
                (0, I03+2000), (cen3-20,cen3+20), (10,120),
                (0, I04+2000), (cen4-30,cen4+30), (10,300),
                (0, I05+2000), (cen5-20,cen5+15), (10,300),
                (0, I0max), (cen6-15, cen6+15), (Gmin,Gmax),
                (m-0.01,m+0.01), (0,np.inf))

        bnds = np.array(bnds)
        bnds = np.transpose(bnds)         
    
        test = dict()
        guesses = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   m, c]
        p_s = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
               'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
               'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
               'm', 'c',]
        
        for i in range(len(p_s)):
            test[p_s[i]] = guesses[i]            
            
            
            
    return(guesses, bnds)

