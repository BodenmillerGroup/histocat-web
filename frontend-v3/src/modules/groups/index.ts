import create from "zustand";
import { schema, normalize } from "normalizr";
import { IGroup, IGroupCreate, IGroupUpdate } from "./models";
import { IMember } from "../members/models";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { isEqual } from "lodash-es";
import { AppToaster } from "../../utils/toaster";

const groupSchema = new schema.Entity("groups");
const groupListSchema = [groupSchema];

type GroupsState = {
  ids: ReadonlyArray<number>;
  entities: { [key: number]: IGroup };
  activeGroupId: number | null;
  myMember: IMember | null;
  groupsTags: string[] | null;

  isGroupAdmin(): boolean;
  getGroups(): Promise<void>;
  getTags(): Promise<void>;
  createGroup(payload: IGroupCreate): Promise<void>;
  getGroup(id: number): Promise<void>;
  updateGroup(id: number, params: IGroupUpdate): Promise<void>;
  deleteGroup(id: number): Promise<void>;
  joinGroup(id: number): Promise<void>;
  getMyMember(groupId: number): Promise<void>;
};

export const useGroupsStore = create<GroupsState>((set, get) => ({
  ids: [],
  entities: {},
  activeGroupId: null,
  myMember: null,
  groupsTags: null,

  isGroupAdmin() {
    const myMember = get().myMember;
    return myMember ? myMember.role >= 100 : false;
  },

  async getGroups() {
    try {
      const data = await api.getGroups();
      if (data) {
        const normalizedData = normalize<IGroup>(data, groupListSchema);
        set({
          ids: normalizedData.result,
          entities: normalizedData.entities.groups ? Object.freeze(normalizedData.entities.groups) : {},
        });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getTags() {
    try {
      const data = await api.getTags();
      if (data && !isEqual(data, get().groupsTags)) {
        set({ groupsTags: data });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async createGroup(params: IGroupCreate) {
    try {
      const data = await api.createGroup(params);
      if (data) {
        set({
          ids: get().ids.concat(data.id),
          entities: Object.freeze({ ...get().entities, [data.id]: data }),
        });
        AppToaster.show({ message: "Group successfully created", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getGroup(id: number) {
    try {
      const data = await api.getGroup(id);
      if (data) {
        const existingId = get().ids.find((id) => id === data.id);
        if (!existingId) {
          set({ ids: get().ids.concat(data.id) });
        }
        set({ entities: Object.freeze({ ...get().entities, [data.id]: data }) });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateGroup(id: number, params: IGroupUpdate) {
    try {
      const data = await api.updateGroup(id, params);
      if (data) {
        set({ entities: Object.freeze({ ...get().entities, [data.id]: data }) });
        AppToaster.show({ message: "Group successfully updated", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async deleteGroup(id: number) {
    try {
      const data = await api.deleteGroup(id);
      if (data) {
        const entities = { ...get().entities };
        delete entities[id];
        set({ ids: get().ids.filter((item) => item !== id), entities: Object.freeze(entities) });
        AppToaster.show({ message: "Group successfully deleted", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async joinGroup(id: number) {
    try {
      const data = await api.joinGroup(id);
      if (data) {
        set({ entities: Object.freeze({ ...get().entities, [data.id]: data }) });
        AppToaster.show({ message: "You successfully joined the group", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getMyMember(groupId: number) {
    try {
      const data = await api.getMyMember(groupId);
      if (data) {
        set({ myMember: data });
      }
    } catch (error) {
      displayApiError(error);
    }
  },
}));
