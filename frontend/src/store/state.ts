import { MainState } from '@/modules/main/state';
import { UserState } from '@/modules/user/state';
import { ExperimentsState } from '@/modules/experiment/state';

export interface State {
  main: MainState;
  user: UserState;
  experiment: ExperimentsState;
}
