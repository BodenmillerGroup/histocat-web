<template>
  <v-container fluid>
    <v-toolbar dense>
      <v-toolbar-title>Edit Project</v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn @click="cancel" text color="primary">Cancel</v-btn>
        <v-btn @click="reset" text color="primary">Reset</v-btn>
        <v-btn @click="submit" text :disabled="!valid" color="primary">Save</v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-card class="mt-4 px-4">
      <v-card-text>
        <v-form v-model="valid" ref="form" lazy-validation>
          <v-text-field label="Name" v-model="name" :rules="nameRules"></v-text-field>
          <v-text-field label="Description" v-model="description"></v-text-field>
          <v-combobox
            v-model="tags"
            :filter="filter"
            :hide-no-data="!search"
            :items="items"
            :search-input.sync="search"
            hide-selected
            label="Tags"
            multiple
            small-chips
            solo
          >
            <template v-slot:no-data>
              <v-list-item>
                <span class="text-subtitle-1">Create</span>
                <v-chip label small>
                  {{ search }}
                </v-chip>
              </v-list-item>
            </template>
            <template v-slot:selection="{ attrs, item, parent, selected }">
              <v-chip v-if="item === Object(item)" v-bind="attrs" :input-value="selected" label small>
                <span class="pr-2">
                  {{ item.text }}
                </span>
                <v-icon small @click="parent.selectItem(item)">mdi-close</v-icon>
              </v-chip>
            </template>
            <template v-slot:item="{ index, item }">
              <v-text-field
                v-if="editing === item"
                v-model="editing.text"
                autofocus
                flat
                background-color="transparent"
                hide-details
                solo
                @keyup.enter="edit(index, item)"
              ></v-text-field>
              <v-chip v-else dark label small>
                {{ item.text }}
              </v-chip>
              <v-spacer></v-spacer>
              <v-list-item-action @click.stop>
                <v-btn icon @click.stop.prevent="edit(index, item)">
                  <v-icon>{{ editing !== item ? "mdi-pencil" : "mdi-check" }}</v-icon>
                </v-btn>
              </v-list-item-action>
            </template>
          </v-combobox>
        </v-form>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { IProjectUpdate } from "@/modules/projects/models";
import { required } from "@/utils/validators";
import { Component, Vue, Watch } from "vue-property-decorator";
import { groupModule } from "@/modules/group";

@Component
export default class EditProject extends Vue {
  readonly groupContext = groupModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  readonly nameRules = [required];

  valid = true;
  name = "";
  description = "";
  tags: any[] = [];

  editing = null;
  index = -1;
  nonce = 1;
  search = null;

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  get items(): any[] {
    const list = this.projectsContext.getters.tags;
    return list.map((item) => {
      return {
        text: item,
      };
    });
  }

  get project() {
    return this.projectsContext.getters.getProject(+this.$router.currentRoute.params.projectId);
  }

  async mounted() {
    await Promise.all([
      this.projectsContext.actions.getProjectsTags(this.activeGroupId!),
      this.projectsContext.actions.getProject(+this.$router.currentRoute.params.projectId),
    ]);
    this.reset();
  }

  reset() {
    this.name = "";
    this.description = "";
    this.tags = [];
    if (this.$refs.form) {
      (this.$refs.form as any).resetValidation();
    }
    if (this.project) {
      this.name = this.project.name;
      this.description = this.project.description;
      this.tags = this.project.tags
        ? this.project.tags.map((item) => {
            return {
              text: item,
            };
          })
        : [];
    }
  }

  cancel() {
    this.$router.back();
  }

  @Watch("tags")
  onTagsChange(val: any[], prev: any[]) {
    if (!val || !prev || val.length === prev.length) {
      return;
    }

    this.tags = val.map((v) => {
      if (typeof v === "string") {
        v = {
          text: v,
        };
        this.items.push(v);
        this.nonce++;
      }
      return v;
    });
  }

  edit(index: number, item) {
    if (!this.editing) {
      this.editing = item;
      this.index = index;
    } else {
      this.editing = null;
      this.index = -1;
    }
  }

  filter(item, queryText: string, itemText: string) {
    const hasValue = (val) => (val != null ? val : "");

    const text = hasValue(itemText);
    const query = hasValue(queryText);

    return text.toString().toLowerCase().indexOf(query.toString().toLowerCase()) > -1;
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const data: IProjectUpdate = {};
      if (this.name) {
        data.name = this.name;
      }
      if (this.description) {
        data.description = this.description;
      }
      if (this.tags && this.tags.length > 0) {
        data.tags = this.tags.map((tag) => tag.text);
      }
      await this.projectsContext.actions.updateProject({ id: this.project!.id, data: data });
      this.$router.back();
    }
  }
}
</script>
