class Point:
    x=0.0		#x coordinate
    y=0.0		#y coordinate
    p=0.0		#potential
    eps=0.0		#dielectric
    cn=0.0		#electon conc
    cp=0.0		#hole conc
    ca=0.0		#acceptor conc
    cd=0.0		#donor conc
    disc_l=0.0		#discretisation for left node
    disc_m=0.0		#discretisation for middle node
    disc_r=0.0		#discretisation for right node
    disc_rhs=0.0	#discretisation for right hand side
    fixed=False		#is this a fixed point

    def fix_point(self,p_fixed):
        self.disc_l=0.0
        self.disc_m=1.0
        self.disc_r=0.0
        self.disc_rhs=p_fixed
        self.p=p_fixed
	self.fixed=True
	print 'Fixing point to',p_fixed
