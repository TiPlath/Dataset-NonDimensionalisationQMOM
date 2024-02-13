function PDF = convertCDFtoPDF(CDF)
    PDF = zeros(length(CDF),2);
    PDF(:,1) = CDF(:,1);
    PDF(:,2) = [diff(CDF(:,2))./diff(CDF(:,1)); 0];