% Convert CDF to PDF using a right Riemann sum.
function CDF = convertPDFtoCDF(PDF)
    CDF = zeros(length(PDF),2);
    CDF(:,1) = PDF(:,1);
    for i = 2:length(PDF)
        CDF(i,2) = CDF(i-1,2) + (PDF(i,1) - PDF(i-1,1)) * PDF(i-1,2);
    end