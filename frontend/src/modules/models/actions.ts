import { IModelUpdate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { ModelsState } from ".";
import { api } from "./api";
import { ModelsGetters } from "./getters";
import { ModelsMutations } from "./mutations";
import { groupModule } from "@/modules/group";

export class ModelsActions extends Actions<ModelsState, ModelsGetters, ModelsMutations, ModelsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
  }

  async getModels(groupId: number) {
    try {
      const data = await api.getGroupModels(groupId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getModel(modelId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const entity = await api.getModel(groupId, modelId);
      if (entity) {
        this.mutations.setEntity(entity);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateModel(payload: { modelId: number; data: IModelUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const entity = await api.updateModel(groupId, payload.modelId, payload.data);
      this.mutations.updateEntity(entity);
      this.main!.mutations.addNotification({ content: "Model successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteModel(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const entity = await api.deleteModel(groupId, id);
      this.mutations.deleteEntity(entity);
      this.main!.mutations.addNotification({ content: "Model successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createModel(formData: FormData) {
    try {
      // await api.uploadValidationFile(
      //   this.main!.getters.token,
      //   payload.validationId,
      //   payload.formData,
      //   () => {
      //     console.log("Upload has started.");
      //     this.main!.mutations.setProcessing(true);
      //   },
      //   () => {
      //     console.log("Upload completed successfully.");
      //     this.main!.mutations.setProcessing(false);
      //     this.main!.mutations.setProcessingProgress(0);
      //     this.main!.mutations.addNotification({ content: "File successfully uploaded", color: "success" });
      //   },
      //   event => {
      //     const percent = Math.round((100 * event.loaded) / event.total);
      //     this.main!.mutations.setProcessingProgress(percent);
      //   },
      //   () => {}
      // );
      const groupId = this.group?.getters.activeGroupId!;
      const entity = await api.createModel(groupId, formData);
      this.mutations.addEntity(entity);
      this.main!.mutations.addNotification({ content: "Model successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
