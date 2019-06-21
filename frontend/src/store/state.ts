import { MainState } from '@/modules/main/state';
import { UserState } from '@/modules/user/state';
import { ExperimentsState } from '@/modules/experiment/state';

export interface RootState {
  main: MainState;
  user: UserState;
  experiment: ExperimentsState;
}
