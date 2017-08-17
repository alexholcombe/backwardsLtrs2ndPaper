function result = pdf_Normal_Exp1(x,p,mu1,sigma1)
    
    pseudo_uniform = [0:1:11 12*ones(1,13) 11:-1:0];
    pseudo_normal = normpdf(-18:1:18,mu1,sigma1).*pseudo_uniform;
    
    normFactor_normal = sum(pseudo_normal);
    
    if normFactor_normal == 0
        normFactor_normal = 0.00001;
    end
    
    uniResultTemp = interp1(-18:1:18,pseudo_uniform,x);
    normResultTemp = normpdf(x,mu1,sigma1).*uniResultTemp;
    
    normResultTemp = normResultTemp/normFactor_normal;
   
    result = p*normResultTemp;

end