function chain = add_joi_to_robot(chain)
%
% Add joints of interest (JOI) with 'box_added' for IK and visualizaiton
%

switch chain.name
    case 'coman'
        chain = add_joi_to_coman(chain);
    case 'iiwa7'
        chain = add_joi_to_iiwa7(chain);
    otherwise
        fprintf(2,'[add_joi_to_robot] unknown [%s].\n',chain.name);
end