C
C.... Table Fill_Ip_Vertex.
#ifdef LONG_POINTER
      Integer*8 Fill_Ip_Vertex
#else
      Integer*4 Fill_Ip_Vertex
#endif
      External  Fill_Ip_Vertex
C
      Integer*4    Fill_Ip_Vertex_ID(1)
      Equivalence( PNT_Common(2),
     +             Fill_Ip_Vertex_ID )
      Integer*4    Fill_Ip_Vertex_flag(1)
      Equivalence( PNT_Common(3),
     +             Fill_Ip_Vertex_flag )
      Real*4       Fill_Ip_Vertex_v(1,3)
      Equivalence( PNT_Common(4),
     +             Fill_Ip_Vertex_v )
      Real*4       Fill_Ip_Vertex_error(1,6)
      Equivalence( PNT_Common(7),
     +             Fill_Ip_Vertex_error )
      Integer*4    Fill_Ip_Vertex_ndf(1)
      Equivalence( PNT_Common(13),
     +             Fill_Ip_Vertex_ndf )
      Real*4       Fill_Ip_Vertex_chisq(1)
      Equivalence( PNT_Common(14),
     +             Fill_Ip_Vertex_chisq )
      Real*4       Fill_Ip_Vertex_cl(1)
      Equivalence( PNT_Common(15),
     +             Fill_Ip_Vertex_cl )
      Real*4       Fill_Ip_Vertex_vx(1)
      Equivalence( Fill_Ip_Vertex_vx,
     +             Fill_Ip_Vertex_v(1,1))
      Real*4       Fill_Ip_Vertex_vy(1)
      Equivalence( Fill_Ip_Vertex_vy,
     +             Fill_Ip_Vertex_v(1,2))
      Real*4       Fill_Ip_Vertex_vz(1)
      Equivalence( Fill_Ip_Vertex_vz,
     +             Fill_Ip_Vertex_v(1,3))
C
