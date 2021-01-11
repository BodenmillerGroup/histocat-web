import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { SegmentationState } from ".";
import { api } from "./api";
import { SegmentationGetters } from "./getters";
import { SegmentationMutations } from "./mutations";
import { groupModule } from "@/modules/group";
import { projectsModule } from "@/modules/projects";
import { ISegmentationSubmission } from "./models";

export class SegmentationActions extends Actions<
  SegmentationState,
  SegmentationGetters,
  SegmentationMutations,
  SegmentationActions
> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  projects?: Context<typeof projectsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.projects = projectsModule.context(store);
  }

  async processSegmentation(params: ISegmentationSubmission) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const projectId = this.projects?.getters.activeProjectId!;
      const data = await api.processSegmentation(groupId, projectId, params);
      this.main!.mutations.addNotification({ content: "Segmentation processing submitted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
