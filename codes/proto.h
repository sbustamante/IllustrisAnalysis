// Input Output Module
int read_head( char * );
int read_snap_gas( char *, int );
int read_snap_all( char *, int );
int ascii_data_all( struct part *, char *, int, int );
int ascii_data_gas( struct part *, char *, int );
int ascii_data_pos( struct part *, char *, int, int );
int ascii_data_slide( struct part *, char *, int, int, int, float, float );

//Gas Module
float pressure( float, float );
float temperature( float );