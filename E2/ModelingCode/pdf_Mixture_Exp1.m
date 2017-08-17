function result = pdf_Mixture_Exp1(x,p,mu1,sigma1)
    
    pseudo_uniform = [0:1:11 12*ones(1,13) 11:-1:0];
    pseudo_normal = normpdf(-18:1:18,mu1,sigma1).*pseudo_uniform;
    
    normFactor_uniform = sum(pseudo_uniform);
    normFactor_normal = sum(pseudo_normal);
    
    if normFactor_uniform == 0
        normFactor_uniform = 0.00001;
    end
    
    if normFactor_normal == 0
        normFactor_normal = 0.00001;
    end
    
    uniResultTemp = interp1(-18:1:18,pseudo_uniform,x);
    normResultTemp = normpdf(x,mu1,sigma1).*uniResultTemp;
    
    uniResultTemp = uniResultTemp/normFactor_uniform;
    normResultTemp = normResultTemp/normFactor_normal;
   
    normResult = p*normResultTemp;
    uniResult = (1-p)*uniResultTemp;
    
    if sum(size(normResult)==size(uniResult))==2
        result = normResult+uniResult;
    else
        result = normResult+uniResult';
    end

end