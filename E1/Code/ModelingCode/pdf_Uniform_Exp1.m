function result = pdf_Uniform_Exp1(x,p)
    
    %fprintf('\np %0.2f, mu1 %2.2f, sigma1 %2.2f', p, mu1, sigma1);
    % Need to taper the normal distribution with the pseudouniform
    
    pseudo_uniform = [0:1:11 12*ones(1,13) 11:-1:0];
    normFactor_uniform = sum(pseudo_uniform);
    
    if normFactor_uniform == 0
        normFactor_uniform = 0.00001;
    end
    
    uniResultTemp = interp1(-18:1:18,pseudo_uniform,x);
    result = p*(uniResultTemp/normFactor_uniform);

end