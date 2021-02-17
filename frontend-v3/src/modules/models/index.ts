import create from "zustand";
import { schema, normalize } from "normalizr";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { AppToaster } from "../../utils/toaster";
import { IModel, IModelUpdate } from "./models";

export const modelSchema = new schema.Entity("models");
export const modelListSchema = [modelSchema];

type ModelsState = {
  ids: ReadonlyArray<number>;
  entities: { [key: number]: IModel };
  activeModelId: number | null;

  setActiveModelId(id: number | null): void;
  getActiveModel(): IModel | null;
  getModels(): Promise<void>;
  getModel(modelId: number): Promise<void>;
  updateModel(modelId: number, params: IModelUpdate): Promise<void>;
  deleteModel(id: number): Promise<void>;
  createModel(formData: FormData): Promise<void>;
};

export const useModelsStore = create<ModelsState>((set, get) => ({
  ids: [],
  entities: {},
  activeModelId: null,

  setActiveModelId(id: number | null) {
    set({ activeModelId: id });
  },

  getActiveModel() {
    const activeModelId = get().activeModelId;
    return activeModelId ? get().entities[activeModelId] : null;
  },

  async getModels() {
    try {
      const data = await api.getAllModels();
      if (data) {
        const normalizedData = normalize<IModel>(data, modelListSchema);
        set({
          ids: normalizedData.result,
          entities: normalizedData.entities.models ? normalizedData.entities.models : {},
        });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getModel(modelId: number) {
    try {
      const data = await api.getModel(modelId);
      if (data) {
        const existingId = get().ids.find((id) => id === data.id);
        if (!existingId) {
          set({ ids: get().ids.concat(data.id) });
        }
        set({ entities: { ...get().entities, [data.id]: data } });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateModel(modelId: number, params: IModelUpdate) {
    try {
      const data = await api.updateModel(modelId, params);
      if (data) {
        set({ entities: { ...get().entities, [data.id]: data } });
        AppToaster.show({ message: "Model successfully updated", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async deleteModel(id: number) {
    try {
      const data = await api.deleteModel(id);
      if (data) {
        const entities = { ...get().entities };
        delete entities[id];
        set({ ids: get().ids.filter((item) => item !== id), entities: entities });
        AppToaster.show({ message: "Model successfully deleted", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

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
      const data = await api.createModel(formData);
      if (data) {
        set({ ids: get().ids.concat(data.id), entities: { ...get().entities, [data.id]: data } });
        AppToaster.show({ message: "Model successfully created", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },
}));
