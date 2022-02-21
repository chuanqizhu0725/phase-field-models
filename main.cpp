#include "domain.h"
#include "fields.h"

int main()
{
    Domain newDomain(1000.0, 1000.0, 1000.0, 1000, 1000, 1000, 1.0, 1.0, 1.0);
    Fields fields(1000, 1000, 1000, 10);
    // newDomain.about();
    cout << fields.PhaseField() << endl;

    return 0;
};