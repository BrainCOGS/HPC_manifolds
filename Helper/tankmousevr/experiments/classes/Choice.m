classdef Choice < uint32
  
  enumeration
    L(1)
    R(2)
    nil(inf)
  end
  
  methods (Static)
    function choices = all()
      choices = enumeration('Choice')';
      choices = choices(1:end-1);
    end
    
    function num = count()
      num = numel(enumeration('Choice'));
    end
  end
  
  methods
    %% UGLY : return the opposite Choice; this logic doesn't work unless there are exactly two non-nil Choice values
    function opp = opposite(obj)
      numValues   = numel(Choice.all());
      assert(numValues == 2);     % the concept of "opposite" only works for sets of 2
      
      flipped     = double(obj);
      flipped     = numValues+1 - flipped;
      opp         = obj;
      sel         = opp >= 1 & opp <= numValues;
      opp(sel)    = flipped(sel);
    end
    
    %% Convert Choice to a nominal sign [-1, 1]
    function sgn = sign(obj)
      sgn         = 2*double(obj) - 3;
      sgn(obj == Choice.nil)  = nan;
    end

  end
  
end
